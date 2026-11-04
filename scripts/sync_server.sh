#!/usr/bin/env bash
# Copy the repository's tracked files to the
# server, using the connection settings of the VS Code SFTP extension (.vscode/sftp.json). Only files that
# differ are sent; nothing on the server is deleted, so server-only files (config/server.yaml,
# config/samples*.tsv, results) are untouched.
# usage: scripts/sync_server.sh [-n] [NAME]
#   -n    show what would be sent
#   NAME  another server, with its settings in .vscode/sftp.NAME.json (same format)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
dry=""; [ "${1:-}" = "-n" ] && { dry="-n"; shift; }
settings=.vscode/sftp.json; [ -n "${1:-}" ] && settings=.vscode/sftp.$1.json
[ -f "$settings" ] || { echo "no $settings" >&2; exit 1; }
read -r host user key remote < <(python3 -c "
import json
c = json.load(open('$settings'))
c = c[0] if isinstance(c, list) else c
print(c['host'], c['username'], c['privateKeyPath'], c['remotePath'])")
files=$(mktemp)
git ls-files > "$files"
rsync $dry -rlci --files-from="$files" -e "ssh -i $key -o BatchMode=yes -o LogLevel=ERROR" ./ "$user@$host:$remote" \
    | grep -v '^\.' || true
rm -f "$files"
