#!/bin/bash
# Build WASM, copy pkg to project root, and optionally commit & push
set -e

echo "=== Building WASM ==="
cd "$(dirname "$0")/core"
wasm-pack build --target web --release

echo "=== Copying pkg to project root ==="
rm -f pkg/.gitignore
cp -r pkg/* ../pkg/

echo "=== Done! pkg/ updated. ==="

if [ "$1" == "push" ]; then
  shift
  msg="${*:-update WASM pkg}"
  cd ..
  git add pkg/ && git commit -m "$msg" && git push
  echo "=== Pushed! ==="
fi
