// Reliability score shown as a colour-coded bar (green = high confidence).
function color(v) {
  if (v >= 0.6) return "#0ca30c";
  if (v >= 0.3) return "#2a78d6";
  return "#eda100";
}

export default function Reliability({ value }) {
  return (
    <div className="rel">
      <div className="track">
        <div className="fill" style={{ width: `${value * 100}%`, background: color(value) }} />
      </div>
      <span className="num">{value.toFixed(2)}</span>
    </div>
  );
}
