import { useState, useEffect } from "react";
import { getHealth, predict } from "./api.js";
import Molecule from "./components/Molecule.jsx";
import TargetsTable from "./components/TargetsTable.jsx";
import SupportingHits from "./components/SupportingHits.jsx";
import SimilarMolecules from "./components/SimilarMolecules.jsx";

const EXAMPLES = [
  { name: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(=O)O" },
  { name: "Ibuprofen", smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O" },
  { name: "Imatinib", smiles: "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C" },
];

// Default number of nearest neighbours (the k slider was removed from the UI;
// k = 10 is the fixed default for now).
const DEFAULT_TOP_K = 10;

function Section({ cls, title, count, children }) {
  return (
    <section className={`card section ${cls}`}>
      <div className="section-head">
        <h2>{title}</h2>
        {count != null && <span className="count">{count}</span>}
      </div>
      {children}
    </section>
  );
}

export default function App() {
  const [health, setHealth] = useState(null);
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [data, setData] = useState(null);
  const [submitted, setSubmitted] = useState("");

  useEffect(() => {
    getHealth().then(setHealth).catch(() => setHealth(null));
  }, []);

  async function run(sm = smiles) {
    const q = sm.trim();
    if (!q) return;
    setLoading(true);
    setError("");
    try {
      const res = await predict(q, DEFAULT_TOP_K);
      setData(res);
      setSubmitted(q);
    } catch (e) {
      setError(e.message);
      setData(null);
    } finally {
      setLoading(false);
    }
  }

  return (
    <>
      <header className="masthead">
        <div className="inner">
          <div className="masthead-top">
            <div>
              <div className="brand">
                <img className="brand-icon" src="/logo.png" alt="MolTarPred logo" />
                <div className="wordmark">MolTarPred</div>
              </div>
              <p>An open-source, early-stage drug-discovery resource for researchers. MolTarPred is a tool for novel target prediction, off-target and polypharmacology exploration, and drug-repurposing leads.</p>
              <p>Enter a molecule (SMILES) and MolTarPred finds its most structurally similar compounds in a curated ChEMBL 36 database, outputting a ranked list of their known protein targets by a reliability score.</p>
            </div>
          </div>

          <div className="masthead-query">
            <div className="query-row">
              <div className="field">
                <div className="field-head">
                  <label htmlFor="smiles">Enter query molecule (SMILES) here:</label>
                  {health && (
                    <div className="badge">
                      <b>{health.n_molecules.toLocaleString()}</b> molecules ·{" "}
                      <b>{health.n_targets.toLocaleString()}</b> targets · MolTarPred v2
                    </div>
                  )}
                </div>
                <input
                  id="smiles"
                  type="text"
                  value={smiles}
                  placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O"
                  onChange={(e) => setSmiles(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && run()}
                />
              </div>
              <button className="btn" onClick={() => run()} disabled={loading || !smiles.trim()}>
                {loading ? <><span className="spinner" /> Predicting…</> : "Predict targets"}
              </button>
            </div>

            <div className="examples">
              <span>Examples:</span>
              {EXAMPLES.map((ex) => (
                <button key={ex.name} className="btn-ghost" onClick={() => { setSmiles(ex.smiles); run(ex.smiles); }}>
                  {ex.name}
                </button>
              ))}
            </div>

            {error && <div className="error">{error}</div>}
          </div>
        </div>
      </header>

      <main className="page">
        {loading && !data && <div className="center-status"><span className="spinner" /> Running MolTarPred…</div>}

        {data && (
          <>
            <Section cls="sec-query" title="Query Molecule">
              <div className="query-preview">
                <div className="mol"><Molecule smiles={submitted} w={220} h={170} alt="query" /></div>
                <div>
                  <div className="subtle" style={{ marginBottom: 4 }}>SMILES</div>
                  <div className="smiles">{submitted}</div>
                </div>
              </div>
            </Section>

            <Section cls="sec-targets" title="Predicted Targets">
              <TargetsTable targets={data.targets} />
            </Section>

            {data.targets.length > 0 && (
              <Section cls="sec-hits" title="Explore Supporting Hits by Target">
                <SupportingHits targets={data.targets} hits={data.top_hits} />
              </Section>
            )}

            <Section cls="sec-similar" title={`Top ${data.top_hits.length} Similar Molecules`}>
              <SimilarMolecules querySmiles={submitted} hits={data.top_hits} />
            </Section>
          </>
        )}

        {!data && !loading && (
          <div className="center-status">Enter a SMILES string or pick an example to get started.</div>
        )}

        <footer>MolTarPred v2 . Powered by ChEMBL</footer>
      </main>
    </>
  );
}
