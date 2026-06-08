# Build WASM and optionally commit & push
param(
  [switch]$Push,
  [string]$Message = "update WASM pkg"
)

Set-Location "$PSScriptRoot\core"
Write-Host "=== Building WASM ===" -ForegroundColor Cyan
wasm-pack build --target web --release --out-dir ..\pkg

Write-Host "=== Done! pkg/ updated. ===" -ForegroundColor Green

if ($Push) {
  Set-Location "$PSScriptRoot"
  git add pkg/
  git commit -m $Message
  git push
  Write-Host "=== Pushed! ===" -ForegroundColor Green
}
