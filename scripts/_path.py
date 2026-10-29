"""讓 scripts/ 在未安裝套件時也能直接執行（把專案根目錄加入 sys.path），並統一以 UTF-8 輸出。"""
import os
import sys

# 輸出被導向檔案或管線時，Windows 會用系統編碼（中文環境為 cp950），印不出 µ、− 等字元
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUTPUTS = os.path.join(ROOT, "outputs")
