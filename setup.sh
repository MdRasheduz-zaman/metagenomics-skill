#!/usr/bin/env bash
set -euo pipefail

echo "🧬 Bootstrapping metagx (metagenomics skill)..."
CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CWD"

# 1. Ensure 'uv' is available
if ! command -v uv &> /dev/null; then
    echo "uv not found. Installing via curl..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    # shellcheck disable=SC1091
    source "$HOME/.local/bin/env" 2>/dev/null || export PATH="$HOME/.local/bin:$PATH"
fi

# 2. Create/sync the virtual environment and install the package (editable) + serve extras
echo "📦 Creating venv and installing metagx..."
uv venv .venv --python 3.11
# shellcheck disable=SC1091
source .venv/bin/activate
uv pip install -e ".[serve]"

# 3. Project directory scaffolding
mkdir -p data results local_databases

# 4. Sanity check the package + registries
python -c "import metagx; print('metagx', metagx.__version__, '— tools:', ', '.join(metagx.list_tools()))"

# 5. Tell the user how to wire MCP into their client
cat <<EOF

✅ Setup complete.

Bioinformatics tools must also be on PATH (install via conda/mamba):
    mamba install -c bioconda kraken2 bracken fastp megahit minimap2 samtools metabat2

Next:
    metagx tools                                   # list pipeline steps
    metagx interview kraken2                        # see what to ask

MCP (Claude Desktop / Cursor) — add to your client config:
    {
      "mcpServers": {
        "metagx": {
          "command": "$CWD/.venv/bin/python",
          "args": ["$CWD/mcp_server.py"]
        }
      }
    }

Web agents (ChatGPT/Gemini/Perplexity):
    $CWD/.venv/bin/uvicorn mcp_server:app --host 0.0.0.0 --port 8000
EOF
