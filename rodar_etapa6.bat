@echo off
REM Etapa 6 - roda calibracao, bateria do SA e gera a tabela comparativa.
REM Pode fechar e rodar de novo: execucoes ja prontas sao puladas.
cd /d "%~dp0"
python -m pip install --quiet scipy matplotlib
python src\download_calibration.py
if not exist results\calibration\best_params.json (
  python src\calibrate_sa.py -j 3 --time-limit 45 --seeds 1-2
)
python src\run_sa.py -j 3 --params results\calibration\best_params.json
python src\make_table_sa.py
echo.
echo PRONTO. Resultados em results\tabela_comparacao.md e results\figs\
pause
