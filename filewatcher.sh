#!/bin/bash
set -euo pipefail

# don't want to add them to dockerfile, hopefully temporary issue that
# can be solved
apt update && apt install -y procps inotify-tools

MEDIA_DIR="/code/media/target_loader_data"
# add the system name to logfile, the command might come from worker
LOG_FILE="$MEDIA_DIR/deletions-$(uname -n).log"

echo "[$(date '+%F %T')] Media watcher started, monitoring $MEDIA_DIR" >> "$LOG_FILE"

mkdir -p "$(dirname "$LOG_FILE")"

inotifywait -m -r -e delete --format '%w%f' "$MEDIA_DIR" | while read -r FILE
do
    {
        echo "[$(date '+%F %T')] File deleted: $FILE"

        # running processes
        echo "--- ps snapshot ---"
        ps -eo pid,ppid,cmd

        # closed, too heavy
        # open files in MEDIA_DIR
        # echo "--- lsof snapshot ---"
        # lsof +D "$MEDIA_DIR" || echo "no open files"

        echo "---------------------------"
    } >> "$LOG_FILE" 2>&1
done
