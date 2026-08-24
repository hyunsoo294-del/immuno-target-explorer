@echo off
chcp 65001 > nul
cd /d "C:\Users\USER\immuno-target-chatbot"
set CI=true
echo Installing required packages if needed...
python -m pip install -r requirements.txt
if errorlevel 1 (
  echo Python is missing or pip failed.
  echo Install Python from https://www.python.org/downloads/
  echo and check "Add python.exe to PATH".
  pause
  exit /b 1
)
echo.
echo 노블타겟 주소: http://127.0.0.1:8502
echo 이 창은 닫지 마세요. 메일 비서는 다른 창에서 8501로 실행합니다.
echo.
start "" "http://127.0.0.1:8502"
python -m streamlit run app.py --server.headless true --server.address 127.0.0.1 --server.port 8502 --browser.gatherUsageStats false
pause
