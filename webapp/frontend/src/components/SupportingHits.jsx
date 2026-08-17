import { useState, useEffect } from "react";
import Molecule from "./Molecule.jsx";

// Lets the user pick a predicted target and see which of the top-k hits support
// it (they are known to be active against it), with structures.
export default function SupportingHits({ targets, hits }) {
  const [selected, setSelected] = useState(targets[0]?.target_chembl_id ?? "");

  useEffect(() => {
    setSelected(targets[0]?.target_chembl_id ?? "");
  }, [targets]);

  if (!targets.length) return null;

  const target = targets.find((t) => t.target_chembl_id === selected) || targets[0];
  const ligandSet = new Set(target.supporting_ligands);
  const supporting = hits.filter((h) => ligandSet.has(h.ligand_chembl_id));

  return (
    <>
      <div className="toolbar">
        <label className="pill" htmlFor="target-select">Target:</label>
        <select id="target-select" value={selected} onChange={(e) => setSelected(e.target.value)}>
          {targets.map((t) => (
            <option key={t.target_chembl_id} value={t.target_chembl_id}>
              {t.target_name} ({t.target_chembl_id})
            </option>
          ))}
        </select>
      </div>

      <p className="subtle" style={{ marginTop: 0 }}>
        <b style={{ color: "var(--text)" }}>{target.target_name}</b> — {supporting.length} of the top hits
        active against this target · reliability {target.reliability.toFixed(2)}
      </p>

      <div className="mol-grid">
        {supporting.map((h, i) => (
          <div className="mol-card" key={h.ligand_chembl_id + "-" + i}>
            <Molecule smiles={h.smiles} w={260} h={200} alt={h.ligand_chembl_id} />
            <div className="cap">
              <b>{h.ligand_chembl_id}</b><br />
              {(h.similarity * 100).toFixed(1)}% similar ·{" "}
              <a className="ext" href={h.chembl_url} target="_blank" rel="noreferrer">ChEMBL ↗</a>
            </div>
          </div>
        ))}
      </div>
    </>
  );
}
