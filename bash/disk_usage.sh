#!/bin/bash

OUTPUT_FILE=".disk_usage"
TEMP_FILE=$(mktemp)
TOTAL_BYTES=0

{   
    echo "Project: $(pwd)"
    echo "Date:    $(date '+%Y-%m-%d %H:%M:%S')"
    echo "================================================================================"
} > "$OUTPUT_FILE"

for item in *; do
    [ -e "$item" ] || continue

    if [ -d "$item" ]; then
        if [ -L "$item" ]; then
            raw_size=$(du -shL "$item" 2>/dev/null | cut -f1)
            size=$(echo "$raw_size" | sed 's/[.,][0-9]//')
            target=$(readlink -f "$item")
            
            echo "${size};[LINK];${item};-> ${target}" >> "$TEMP_FILE"
            bytes=$(du -sbL "$item" 2>/dev/null | cut -f1)
        else
            raw_size=$(du -sh "$item" 2>/dev/null | cut -f1)
            size=$(echo "$raw_size" | sed 's/[.,][0-9]//')
            
            echo "${size};[DIR];${item}; " >> "$TEMP_FILE"
            bytes=$(du -sb "$item" 2>/dev/null | cut -f1)
        fi

        if [[ "$bytes" =~ ^[0-9]+$ ]]; then
            TOTAL_BYTES=$((TOTAL_BYTES + bytes))
        fi
    fi
done

HUMAN_TOTAL=$(numfmt --to-unit=1G --format="%.0fG" $TOTAL_BYTES 2>/dev/null || echo "0G")

{
    echo -e "Size;Type;Location;" | column -s ";" -t
    echo "--------------------------------------------------------------------------------"
    column -s ";" -t "$TEMP_FILE"
    echo "--------------------------------------------------------------------------------"
    echo -e "${HUMAN_TOTAL};[TOTAL]; ; " | column -s ";" -t
    echo "================================================================================"
} >> "$OUTPUT_FILE"

rm "$TEMP_FILE"
cat "$OUTPUT_FILE"
