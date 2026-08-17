import { moleculeUrl } from "../api.js";

// Renders a molecule as a backend-generated SVG. Transparent background blends
// with whatever card it sits in.
export default function Molecule({ smiles, w = 300, h = 220, alt = "molecule" }) {
  return (
    <img
      src={moleculeUrl(smiles, w, h)}
      width={w}
      height={h}
      alt={alt}
      loading="lazy"
      style={{ maxWidth: "100%", height: "auto" }}
    />
  );
}
