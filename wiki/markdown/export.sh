#/bin/bash

pushd "$(dirname "$0")"

for f in *.md; do
  name=${f%.md}
  pandoc --template=template.html \
    --katex \
    --css=wiki.css \
    --lua-filter=urlfilter.lua \
    --metadata title="$name" \
    -i "$name.md" \
    -o "$name.html"
  mv "$name.html" ../export/
done

popd
