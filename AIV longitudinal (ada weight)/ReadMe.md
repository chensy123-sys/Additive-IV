```bash

n=5000
nohup Rscript script/run_simulation.R $n 0,0 > logs/mt$n_00.log 2>&1 &
nohup Rscript script/run_simulation.R $n 0,1 > logs/mt$n_01.log 2>&1 &
nohup Rscript script/run_simulation.R $n 1,0 > logs/mt$n_10.log 2>&1 &
nohup Rscript script/run_simulation.R $n 1,1 > logs/mt$n_11.log 2>&1 &

nohup Rscript script/run_simulation.R 5000 0,0 > logs/mt_5000_00.log 2>&1 &
nohup Rscript script/run_simulation.R 5000 1,0 > logs/mt_5000_10.log 2>&1 &
nohup Rscript script/run_simulation.R 5000 0,1 > logs/mt_5000_01.log 2>&1 &
nohup Rscript script/run_simulation.R 5000 1,1 > logs/mt_5000_11.log 2>&1 &
```

ps -ef | grep run_simulation