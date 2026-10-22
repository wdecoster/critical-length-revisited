#!/usr/bin/env bash
# Copy the repository's tracked files to the
# server, using the connection settings of the VS Code SFTP extension (.vscode/sftp.json). Only files that
# differ are sent; nothing on the server is deleted, so server-only files (config/server.yaml,
# config/samples*.tsv, results) are untouched. -n: show what would be sent.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
dry=""; [ "${1:-}" = "-n" ] && dry="-n"
read -r host user key remote < <(python3 -c "
import json
c = json.load(open('.vscode/sftp.json'))
c = c[0] if isinstance(c, list) else c
print(c['host'], c['username'], c['privateKeyPath'], c['remotePath'])")
files=$(mktemp)
git ls-files > "$files"
rsync $dry -rlci --files-from="$files" -e "ssh -i $key -o BatchMode=yes -o LogLevel=ERROR" ./ "$user@$host:$remote" \
    | grep -v '^\.' || true
rm -f "$files"
