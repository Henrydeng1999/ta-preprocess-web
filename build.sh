#!/bin/bash
# Build WASM and optionally commit & push
set -e

echo "=== Building WASM ==="
cd "$(dirname "$0")/core"
wasm-pack build --target web --release --out-dir ../pkg

echo "=== Done! pkg/ updated. ==="

if [ "$1" == "push" ]; then
  shift
  msg="${*:-update WASM pkg}"
  cd ..
  git add pkg/ && git commit -m "$msg" && git push
  echo "=== Pushed! ==="
fi
