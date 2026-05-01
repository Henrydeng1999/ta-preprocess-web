use crate::TaData;

fn read_be_uint32(data: &[u8], offset: usize) -> Option<u32> {
    if offset + 4 > data.len() { return None; }
    Some(u32::from_be_bytes([data[offset], data[offset + 1], data[offset + 2], data[offset + 3]]))
}

fn read_be_uint16(data: &[u8], offset: usize) -> Option<u16> {
    if offset + 2 > data.len() { return None; }
    Some(u16::from_be_bytes([data[offset], data[offset + 1]]))
}

fn read_be_float64(data: &[u8], offset: usize) -> Option<f64> {
    if offset + 8 > data.len() { return None; }
    Some(f64::from_be_bytes([
        data[offset], data[offset + 1], data[offset + 2], data[offset + 3],
        data[offset + 4], data[offset + 5], data[offset + 6], data[offset + 7],
    ]))
}

fn read_ufs_string(data: &[u8], offset: usize, length: usize) -> String {
    let end = std::cmp::min(offset + length, data.len());
    let mut result = String::new();
    for i in offset..end {
        let b = data[i];
        if b >= 32 && b <= 126 {
            result.push(b as char);
        } else {
            break;
        }
    }
    result
}

struct SectionHeader {
    count: usize,
    data_offset: usize,
}

fn parse_section_header(data: &[u8], mut offset: usize) -> Result<SectionHeader, String> {
    let name_len = read_be_uint32(data, offset).ok_or("UFS: failed to read section name length")? as usize;
    offset += 4;
    if offset + name_len > data.len() {
        return Err("UFS: section name exceeds file".to_string());
    }
    offset += name_len;
    offset += 4; // pad + dtype
    offset += 2; // unit
    let count = read_be_uint32(data, offset).ok_or("UFS: failed to read section count")? as usize;
    offset += 4;
    Ok(SectionHeader { count, data_offset: offset })
}

struct IntensityHeader {
    wl_count: usize,
    t_count: usize,
    data_offset: usize,
}

fn parse_intensity_header(data: &[u8], mut offset: usize) -> Result<IntensityHeader, String> {
    let name_len = read_be_uint32(data, offset).ok_or("UFS: failed to read intensity name length")? as usize;
    offset += 4;
    if offset + name_len > data.len() {
        return Err("UFS: intensity name exceeds file".to_string());
    }
    offset += name_len;
    offset += 6; // pad + dtype + unit
    let wl_count = read_be_uint16(data, offset).ok_or("UFS: failed to read wl_count")? as usize;
    offset += 2;
    offset += 2; // pad
    let t_count = read_be_uint16(data, offset).ok_or("UFS: failed to read t_count")? as usize;
    offset += 2;
    Ok(IntensityHeader { wl_count, t_count, data_offset: offset })
}

pub fn parse_ufs_file(data: &[u8]) -> Result<TaData, String> {
    if data.len() < 16 {
        return Err("File too small".to_string());
    }

    let mut offset = 0;
    let ver_len = read_be_uint32(data, offset).ok_or("UFS: failed to read version length")? as usize;
    offset += 4;
    if offset + ver_len > data.len() {
        return Err("UFS: version string exceeds file".to_string());
    }
    let version = read_ufs_string(data, offset, ver_len);
    if version != "Version2" {
        return Err(format!("Not a valid UFS file (got: {})", version));
    }
    offset += ver_len;

    let wl_section = parse_section_header(data, offset)
        .map_err(|e| format!("UFS wavelength section: {}", e))?;
    let mut wavelengths = Vec::with_capacity(wl_section.count);
    for i in 0..wl_section.count {
        let off = wl_section.data_offset + i * 8;
        match read_be_float64(data, off) {
            Some(v) => wavelengths.push(v),
            None => return Err(format!("UFS: wavelength[{}] out of bounds", i)),
        }
    }
    offset = wl_section.data_offset + wl_section.count * 8;

    let t_section = parse_section_header(data, offset)
        .map_err(|e| format!("UFS time section: {}", e))?;
    let mut times = Vec::with_capacity(t_section.count);
    for i in 0..t_section.count {
        let off = t_section.data_offset + i * 8;
        match read_be_float64(data, off) {
            Some(v) => times.push(v),
            None => return Err(format!("UFS: time[{}] out of bounds", i)),
        }
    }
    offset = t_section.data_offset + t_section.count * 8;

    let int_section = parse_intensity_header(data, offset)
        .map_err(|e| format!("UFS intensity section: {}", e))?;
    let wl_count = int_section.wl_count;
    let t_count = int_section.t_count;

    let mut intensity: Vec<Vec<f64>> = Vec::with_capacity(wl_count);
    for wi in 0..wl_count {
        let mut row = Vec::with_capacity(t_count);
        for ti in 0..t_count {
            let flat_idx = wi * t_count + ti;
            let off = int_section.data_offset + flat_idx * 8;
            match read_be_float64(data, off) {
                Some(val) if val.is_finite() && val.abs() >= 1e-300 => row.push(val),
                _ => row.push(f64::NAN),
            }
        }
        intensity.push(row);
    }

    Ok(TaData {
        time_array: times,
        wavelength_array: wavelengths,
        ta_2d: intensity,
    })
}
