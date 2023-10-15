#/bin/bash

pushd "$(dirname "$0")"

for f in *.md; do
  name=${f%.md}
  if [[ $name == "Landing" ]]
  then
    toc=""
  else
    toc="--toc"
  fi
  pandoc \
    --template=../templates/template.html \
    --shift-heading-level-by=1 \
    $toc \
    --katex \
    --css=wiki.css \
    --lua-filter=urlfilter.lua \
    --from markdown+yaml_metadata_block \
    -i "$name.md" \
    -o "$name.html"
  mv "$name.html" ../export/
done

popd
