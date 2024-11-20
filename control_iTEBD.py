import subprocess
import numpy as np
import json

def run_iTEBD(Jz, Jxy, hx, m, D, chi_max, tau_max, tau_min, tau_step, inv_precision, use_complex, second_ST, output_dyn, output_dyn_num):
    # コマンドを構築
    command = [
        "python3", "iTEBD.py",
        "-Jz", str(Jz),
        "-Jxy", str(Jxy),
        "-hx", str(hx),
        "-m", str(m),
        "-D", str(D),
        "-chi", str(chi_max),
        "-tau_max", str(tau_max),
        "-tau_min", str(tau_min),
        "-tau_step", str(tau_step),
        "-inv_precision", str(inv_precision),
    ]
    
    if use_complex:
        command.append("--use_complex")
    if second_ST:
        command.append("--second_ST")
    if output_dyn:
        command.append("--output_dyn")
        command.extend(["-output_dyn_num", str(output_dyn_num)])

    # サブプロセスで実行し、結果を取得
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return result.stdout.strip()  # 標準出力を返す

if __name__ == "__main__":
    # 固定パラメータ
    Jz = 2.1
    Jxy = 1.0
    m = 2
    D = 0.0
    chi_max = 10
    tau_max = 0.1
    tau_min = 0.001
    tau_step = 2000
    inv_precision = 1e-10
    use_complex = True
    second_ST = False
    output_dyn = True
    output_dyn_num = 100

    # hx を 0 から 5 まで 0.1 刻みで変化させる
    hx_values = np.arange(0, 5.1, 0.1)
    results = []

    for hx in hx_values:
        print(f"Running iTEBD simulation for hx={hx:.1f}...")
        output = run_iTEBD(Jz, Jxy, hx, m, D, chi_max, tau_max, tau_min, tau_step, inv_precision, use_complex, second_ST, output_dyn, output_dyn_num)
        results.append((hx, output))
    final_result = {}
    # 結果を表示
    print("\nResults:")
    for hx, result in results:
        print(f"hx={hx:.1f} -> {result}")
        final_result[hx] = result
    with open("results_iTEBD_chi10_without_sy.json", "w") as json_file:
        json.dump(final_result, json_file, indent=4)  # 見やすくインデントを追加
