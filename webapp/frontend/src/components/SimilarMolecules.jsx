import Molecule from "./Molecule.jsx";
import { downloadCsv } from "../api.js";

export default function SimilarMolecules({ querySmiles, hits }) {
  if (!hits.length) return null;

  return (
    <>
      <div className="toolbar">
        <span className="pill">Nearest neighbours by Morgan/Tanimoto similarity</span>
        <button
          className="btn-ghost"
          style={{ marginLeft: "auto" }}
          onClick={() =>
            downloadCsv("moltarpred_similar_molecules.csv", hits, [
              { label: "Ligand ChEMBL ID", get: (h) => h.ligand_chembl_id },
              { label: "SMILES", get: (h) => h.smiles },
              { label: "Tanimoto similarity", get: (h) => h.similarity },
            ])
          }
        >
          ↓ Download CSV
        </button>
      </div>

      <div className="table-wrap" style={{ marginBottom: 18 }}>
        <table id="similar-table">
          <thead>
            <tr>
              <th>#</th>
              <th>Ligand</th>
              <th>SMILES</th>
              <th>Tanimoto</th>
              <th>ChEMBL</th>
            </tr>
          </thead>
          <tbody>
            {hits.map((h, i) => (
              <tr key={h.ligand_chembl_id + i}>
                <td className="subtle">{i + 1}</td>
                <td className="mono">{h.ligand_chembl_id}</td>
                <td className="mono subtle" style={{ maxWidth: 320, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{h.smiles}</td>
                <td className="mono">{h.similarity.toFixed(4)}</td>
                <td><a className="ext" href={h.chembl_url} target="_blank" rel="noreferrer">View ↗</a></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="mol-grid">
        <div className="mol-card query">
          <Molecule smiles={querySmiles} w={260} h={200} alt="query molecule" />
          <div className="cap"><b>Query molecule</b></div>
        </div>
        {hits.map((h, i) => (
          <div className="mol-card" key={h.ligand_chembl_id + "-" + i}>
            <Molecule smiles={h.smiles} w={260} h={200} alt={h.ligand_chembl_id} />
            <div className="cap"><b>{h.ligand_chembl_id}</b><br />{(h.similarity * 100).toFixed(1)}% similar</div>
          </div>
        ))}
      </div>
    </>
  );
}
