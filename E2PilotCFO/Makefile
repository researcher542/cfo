NOW=`date +%Y.%m.%d.%H:%M:%S`
HOSTNAME=`hostname`
main:
	nohup julia --project=. -t auto paper/main.jl &

cs_net:
	nohup julia --project=. -t 128 paper/get_cs_net.jl > log/info.cs_net.$(NOW).$(HOSTNAME).log &

# Get all results for one src-des pair
ipair=364
region=oldmap ## test4, oldmap, TXOH
getone:
	julia --project=. -t 8 paper/get_one_src_des.jl --idx $(ipair) --region $(region) > log/info.getone.$(NOW).$(HOSTNAME).log 

run_nohup:
	nohup python3 run_on_server.py --nthread 8 --startidx 201 --ninstance 210 --region largedata --mem 512 > log/info.main.$(NOW).$(HOSTNAME).log &


# python3 run_on_server.py --tratio --nthread 8 --startidx 0 --ninstance 100 --region largedata > log/info.main.$(NOW).$(HOSTNAME).log
run:
	python3 run_on_server.py --tratio --alpha --future --noise --nthread 8 --startidx 301 --ninstance 310 --region largedata > log/info.main.$(NOW).$(HOSTNAME).log

run_cmp:
	julia --project=. -t 8 paper/method_cmp.jl --dsc true --tbexnet false --grb true 

test:
	echo ${HOSTNAME}

	