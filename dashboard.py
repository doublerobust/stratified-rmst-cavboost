#!/usr/bin/env python3
"""Stratified RMST Simulation Dashboard — real-time progress monitor."""
import csv, glob, json, os, sys, time
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse

REPO = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(REPO, "results_csv")
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8081
SCENARIOS = ["S1_Linear", "S2_Diff", "S3_U", "S4_Enclave", "S5_S", "S6_Cross"]
N_REPS = 50
START_TIME = time.time()
REP_TIMES = []  # rolling window of rep completion times


def get_progress():
    """Count CSV files per scenario. Only counts Sc{N}_Rep{N}.csv format."""
    counts = {i: 0 for i in range(1, 7)}
    total = 0
    try:
        for f in sorted(os.listdir(CSV_DIR)):
            if f.startswith("Sc") and "_Rep" in f and f.endswith(".csv"):
                parts = f.split("_")
                sc = int(parts[0].replace("Sc", ""))
                if 1 <= sc <= 6:
                    counts[sc] += 1
                    total += 1
    except FileNotFoundError:
        pass
    return counts, total


def get_summary():
    """Aggregate results per scenario."""
    results = {}
    try:
        for f in glob.glob(os.path.join(CSV_DIR, "*.csv")):
            with open(f) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    sc = int(row.get("scenario", row.get("sc", 0)))
                    if sc not in results:
                        results[sc] = {"orig": [], "strat": [], "vt": []}
                    results[sc]["orig"].append(float(row["orig_auc"]))
                    results[sc]["strat"].append(float(row["strat_auc"]))
                    results[sc]["vt"].append(float(row["vt_auc"]))
    except (FileNotFoundError, StopIteration):
        pass
    summary = {}
    for sc, vals in results.items():
        n = len(vals["orig"])
        if n == 0:
            continue
        summary[SCENARIOS[sc-1]] = {
            "n": n,
            "orig_auc": sum(vals["orig"]) / n,
            "strat_auc": sum(vals["strat"]) / n,
            "vt_auc": sum(vals["vt"]) / n,
            "orig_correct": sum(1 for v in vals["orig"] if v > 0.5) / n * 100,
            "strat_correct": sum(1 for v in vals["strat"] if v > 0.5) / n * 100,
            "vt_correct": sum(1 for v in vals["vt"] if v > 0.5) / n * 100,
        }
    return summary


def compute_eta(total):
    """Simple ETA based on recent completion rate."""
    elapsed = time.time() - START_TIME
    if total == 0:
        return "calculating..."
    rate = total / elapsed  # reps per second
    if rate == 0:
        return "calculating..."
    remaining = (300 - total) / rate
    if remaining < 60:
        return f"{remaining:.0f}s"
    elif remaining < 3600:
        return f"{remaining/60:.0f}m"
    else:
        return f"{remaining/3600:.1f}h"


HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RMST Simulation Dashboard</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;background:#0d1117;color:#c9d1d9;padding:20px;min-height:100vh}
h1{font-size:1.5rem;margin-bottom:4px;color:#f0f6fc}
.subtitle{color:#8b949e;font-size:0.9rem;margin-bottom:20px}
.status-bar{display:flex;gap:16px;flex-wrap:wrap;margin-bottom:24px}
.stat-card{background:#161b22;border:1px solid #30363d;border-radius:8px;padding:16px 20px;flex:1;min-width:140px}
.stat-card .label{font-size:0.8rem;color:#8b949e;text-transform:uppercase}
.stat-card .value{font-size:1.6rem;font-weight:600;margin-top:4px}
.scenario-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:12px}
.scenario-card{background:#161b22;border:1px solid #30363d;border-radius:8px;padding:16px}
.scenario-card h3{font-size:1rem;margin-bottom:8px;color:#f0f6fc}
.bar-bg{background:#21262d;border-radius:4px;height:20px;overflow:hidden;margin-bottom:4px}
.bar-fill{height:100%;border-radius:4px;transition:width .5s ease;background:linear-gradient(90deg,#238636,#2ea043)}
.bar-fill.done{background:linear-gradient(90deg,#1f6feb,#58a6ff)}
.bar-label{display:flex;justify-content:space-between;font-size:0.8rem;color:#8b949e;margin-bottom:12px}
table{width:100%;border-collapse:collapse;margin-top:16px;font-size:0.85rem}
th{text-align:left;padding:8px 12px;border-bottom:2px solid #30363d;color:#8b949e;font-size:0.75rem;text-transform:uppercase}
td{padding:8px 12px;border-bottom:1px solid #21262d}
.positive{color:#3fb950}
.negative{color:#f85149}
.elapsed{color:#8b949e;font-size:0.75rem;margin-top:12px}
.hidden{display:none}
</style>
</head>
<body>
<h1>🕷️ Stratified RMST Simulation</h1>
<p class="subtitle">50 reps × 6 scenarios — parallel on agent server</p>

<div class="status-bar">
  <div class="stat-card"><div class="label">Progress</div><div class="value" id="total-progress">0/300</div></div>
  <div class="stat-card"><div class="label">Percent</div><div class="value" id="percent">0%</div></div>
  <div class="stat-card"><div class="label">ETA</div><div class="value" id="eta">—</div></div>
  <div class="stat-card"><div class="label">Elapsed</div><div class="value" id="elapsed">0m</div></div>
</div>

<div class="scenario-grid" id="scenarios"></div>

<div id="summary-section" class="hidden">
  <h2 style="margin:24px 0 12px;font-size:1.1rem">📊 Results Summary</h2>
  <table><thead>
    <tr><th>Scenario</th><th>N</th><th>Orig AUC</th><th>Strat AUC</th><th>VT AUC</th><th>Gain</th><th>% Correct (Orig)</th><th>% Correct (Strat)</th></tr>
  </thead><tbody id="summary-body"></tbody></table>
</div>

<p class="elapsed" id="footer"></p>

<script>
const SCENARIOS = ["S1_Linear","S2_Diff","S3_U","S4_Enclave","S5_S","S6_Cross"];
const COLORS = ["#238636","#1f6feb","#d29922","#f85149","#bc8cff","#f0883e"];
const START = Date.now();

function pad(n){return n<10?"0"+n:n}
function fmtTime(ms){let s=Math.floor(ms/1000);let m=Math.floor(s/60);s=s%60;return m+"m "+s+"s"}

function renderBars(counts,total){
  let html="";
  for(let i=0;i<6;i++){
    let n=counts[i+1]||0,pct=n/50*100;
    let cl=pct>=100?"bar-fill done":"bar-fill";
    html+='<div class="scenario-card"><h3>'+SCENARIOS[i]+'</h3>'+
      '<div class="bar-label"><span>'+n+'/50</span><span>'+pct.toFixed(0)+'%</span></div>'+
      '<div class="bar-bg"><div class="'+cl+'" style="width:'+pct+'%"></div></div></div>';
  }
  document.getElementById("scenarios").innerHTML=html;
}

function renderSummary(data){
  let html="";
  for(let sc of SCENARIOS){
    if(!data[sc]){html+="<tr><td>"+sc+"</td><td colspan='7'>running...</td></tr>";continue;}
    let d=data[sc],gain=d.strat_auc-d.orig_auc;
    let gc=d.gain||0;
    html+="<tr><td>"+sc+"</td><td>"+d.n+"</td>"+
      "<td>"+d.orig_auc.toFixed(4)+"</td>"+
      "<td class='positive'>"+d.strat_auc.toFixed(4)+"</td>"+
      "<td>"+d.vt_auc.toFixed(4)+"</td>"+
      "<td class='"+(gain>=0?"positive":"negative")+"'>"+(gain>=0?"+":"")+gain.toFixed(4)+"</td>"+
      "<td>"+d.orig_correct.toFixed(0)+"%</td>"+
      "<td>"+d.strat_correct.toFixed(0)+"%</td></tr>";
  }
  document.getElementById("summary-body").innerHTML=html;
}

async function poll(){
  try{
    let r=await fetch("/api/progress");
    let d=await r.json();
    let total=d.total;
    document.getElementById("total-progress").textContent=total+"/300";
    document.getElementById("percent").textContent=(total/300*100).toFixed(0)+"%";
    document.getElementById("eta").textContent=d.eta;
    document.getElementById("elapsed").textContent=fmtTime(Date.now()-START);
    renderBars(d.counts,total);

    if(total>0){
      let r2=await fetch("/api/summary");
      let s=await r2.json();
      renderSummary(s);
      if(total>=300) document.getElementById("summary-section").classList.remove("hidden");
    }
  }catch(e){console.error("poll error",e)}
  setTimeout(poll,3000);
}
poll();
</script>
</body>
</html>"""


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/progress":
            counts, total = get_progress()
            self._json({"counts": counts, "total": total, "eta": compute_eta(total)})
        elif path == "/api/summary":
            self._json(get_summary())
        else:
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            self.wfile.write(HTML.encode())

    def _json(self, data):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(data).encode())

    def log_message(self, *a):
        pass  # quiet


if __name__ == "__main__":
    print(f"🚀 Dashboard at http://0.0.0.0:{PORT}")
    print(f"   Progress API: http://localhost:{PORT}/api/progress")
    print(f"   Summary API: http://localhost:{PORT}/api/summary")
    HTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
