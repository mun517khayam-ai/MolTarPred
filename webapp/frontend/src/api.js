// Thin API client for the MolTarPred v2 backend.

export async function getHealth() {
  const r = await fetch("/api/health");
  if (!r.ok) throw new Error("Backend unavailable");
  return r.json();
}

export async function predict(smiles, topK) {
  const r = await fetch("/api/predict", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles, top_k: topK }),
  });
  if (r.status === 422) throw new Error("Invalid SMILES string — please check your input.");
  if (!r.ok) throw new Error("Prediction failed. Is the backend running?");
  return r.json();
}

// Returns a same-origin URL for a rendered molecule SVG.
export function moleculeUrl(smiles, w = 320, h = 240, legend = "") {
  const p = new URLSearchParams({ smiles, w, h });
  if (legend) p.set("legend", legend);
  return `/api/molecule?${p.toString()}`;
}

// Build a CSV string + trigger a browser download.
export function downloadCsv(filename, rows, columns) {
  const header = columns.map((c) => c.label).join(",");
  const body = rows
    .map((row) =>
      columns
        .map((c) => {
          const v = c.get(row) ?? "";
          const s = String(v).replace(/"/g, '""');
          return /[",\n]/.test(s) ? `"${s}"` : s;
        })
        .join(",")
    )
    .join("\n");
  const blob = new Blob([header + "\n" + body], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}
