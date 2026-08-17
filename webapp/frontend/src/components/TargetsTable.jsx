import { useState, useEffect } from "react";
import Reliability from "./Reliability.jsx";
import { downloadCsv } from "../api.js";

const PAGE = 10;

// Click a column header to sort. Reliability (descending) is the default.
// Results are paginated 10 per page.
export default function TargetsTable({ targets }) {
  const [sortKey, setSortKey] = useState("rel");
  const [sortDir, setSortDir] = useState(-1); // 1 = asc, -1 = desc
  const [page, setPage] = useState(0);

  // New prediction → back to the first page.
  useEffect(() => {
    setPage(0);
  }, [targets]);

  if (!targets.length) return <p className="subtle">No targets found for this molecule.</p>;

  function onSort(key) {
    if (key === sortKey) setSortDir((d) => -d);
    else {
      setSortKey(key);
      setSortDir(key === "rel" ? -1 : 1);
    }
    setPage(0);
  }

  const value = (t, key) =>
    key === "rel"
      ? t.reliability
      : key === "name"
      ? (t.target_name || "").toLowerCase()
      : (t.organism || "").toLowerCase();

  const sorted = [...targets].sort((a, b) => {
    const va = value(a, sortKey), vb = value(b, sortKey);
    if (va < vb) return -1 * sortDir;
    if (va > vb) return 1 * sortDir;
    return 0;
  });

  const pages = Math.ceil(sorted.length / PAGE);
  const cur = Math.min(page, pages - 1);
  const start = cur * PAGE;
  const end = Math.min(start + PAGE, sorted.length);
  const pageRows = sorted.slice(start, end);

  const arrow = (key) => (key === sortKey ? (sortDir === 1 ? "▲" : "▼") : "↕");
  const th = (key) => ({
    className: `sortable${key === sortKey ? " active" : ""}`,
    onClick: () => onSort(key),
    "aria-sort": key === sortKey ? (sortDir === 1 ? "ascending" : "descending") : "none",
  });

  return (
    <>
      <div className="toolbar">
        <span className="pill">Ranked by reliability score (fraction of top-k hits sharing the target)</span>
        <button
          className="btn-ghost"
          style={{ marginLeft: "auto" }}
          onClick={() =>
            downloadCsv("moltarpred_targets.csv", sorted, [
              { label: "Target ChEMBL ID", get: (t) => t.target_chembl_id },
              { label: "Target name", get: (t) => t.target_name },
              { label: "Organism", get: (t) => t.organism },
              { label: "Reliability", get: (t) => t.reliability },
              { label: "Supporting hits", get: (t) => t.supporting_hits },
            ])
          }
        >
          ↓ Download CSV
        </button>
      </div>
      <div className="table-wrap">
        <table id="targets-table">
          <thead>
            <tr>
              <th>#</th>
              <th {...th("name")}>Target<span className="arrow">{arrow("name")}</span></th>
              <th {...th("org")}>Organism<span className="arrow">{arrow("org")}</span></th>
              <th {...th("rel")}>Reliability<span className="arrow">{arrow("rel")}</span></th>
              <th>ChEMBL</th>
            </tr>
          </thead>
          <tbody>
            {pageRows.map((t, i) => (
              <tr key={t.target_chembl_id}>
                <td className="subtle">{start + i + 1}</td>
                <td>
                  <div>{t.target_name}</div>
                  <div className="mono subtle">{t.target_chembl_id}</div>
                </td>
                <td>{t.organism}</td>
                <td><Reliability value={t.reliability} /></td>
                <td><a className="ext" href={t.chembl_url} target="_blank" rel="noreferrer">View ↗</a></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {pages > 1 && (
        <div className="pager">
          <span className="pager-showing">Showing {start + 1}–{end} of {sorted.length}</span>
          <button
            className="pager-btn"
            aria-label="Previous page"
            disabled={cur === 0}
            onClick={() => setPage((p) => Math.max(0, p - 1))}
          >
            &lsaquo;
          </button>
          <span className="pager-info">{cur + 1} / {pages}</span>
          <button
            className="pager-btn"
            aria-label="Next page"
            disabled={cur === pages - 1}
            onClick={() => setPage((p) => Math.min(pages - 1, p + 1))}
          >
            &rsaquo;
          </button>
        </div>
      )}
    </>
  );
}
