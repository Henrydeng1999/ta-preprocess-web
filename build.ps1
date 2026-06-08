# Build WASM, copy pkg to project root, and optionally commit & push
param(
  [switch]$Push,
  [string]$Message = "update WASM pkg"
)

Set-Location "$PSScriptRoot\core"
Write-Host "=== Building WASM ===" -ForegroundColor Cyan
wasm-pack build --target web --release

Write-Host "=== Copying pkg to project root ===" -ForegroundColor Cyan
Remove-Item "$PSScriptRoot\core\pkg\.gitignore" -ErrorAction SilentlyContinue
Copy-Item "$PSScriptRoot\core\pkg\*" "$PSScriptRoot\pkg\" -Force

Write-Host "=== Done! pkg/ updated. ===" -ForegroundColor Green

if ($Push) {
  Set-Location "$PSScriptRoot"
  git add pkg/
  git commit -m $Message
  git push
  Write-Host "=== Pushed! ===" -ForegroundColor Green
}
