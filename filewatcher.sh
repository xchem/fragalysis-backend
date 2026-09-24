#!/bin/bash
set -euo pipefail

# don't want to add them to dockerfile, hopefully temporary issue that
# can be solved
apt update && apt install -y procps inotify-tools

MEDIA_DIR="/code/media/target_loader_data"
# add the system name to logfile, the command might come from worker
LOG_FILE="$MEDIA_DIR/deletions-$(uname -n).log"

# Deletions we asked for - deleting a target, or peeling an upload off one -
# drop a marker here for as long as they run, and are not logged: this watcher
# is here to catch files vanishing when nothing should have touched them, and a
# target deletion would otherwise bury that in thousands of lines. Written by
# viewer.media_watcher; they sit beside the watched directory rather than inside
# it, so creating and removing them raises no events of their own.
PAUSE_GLOB="$(dirname "$MEDIA_DIR")/.filewatcher-pause."*

deletion_expected() {
    compgen -G "$PAUSE_GLOB" > /dev/null
}

echo "[$(date '+%F %T')] Media watcher started, monitoring $MEDIA_DIR" >> "$LOG_FILE"

mkdir -p "$(dirname "$LOG_FILE")"

suppressed=0

inotifywait -m -r -e delete --format '%w%f' "$MEDIA_DIR" | while read -r FILE
do
    if deletion_expected; then
        suppressed=$((suppressed + 1))
        continue
    fi

    {
        # Say how much was hidden, so a deliberate deletion still leaves a trace
        # of its size without listing every file.
        if [ "$suppressed" -gt 0 ]; then
            echo "[$(date '+%F %T')] ($suppressed expected deletion(s) not logged)"
            suppressed=0
        fi
        echo "[$(date '+%F %T')] File deleted: $FILE"
    } >> "$LOG_FILE" 2>&1
done
