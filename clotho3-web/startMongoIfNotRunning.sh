# this script checks if the mongod is running, starts it if not

if pgrep -q mongod; then
	echo running;
else
	mongod;
fi

exit 0;

#OLD VERSION BELOW

`ps -A | grep -q '[m]ongod'`

if [ "$?" -eq "0" ]; then
    echo "running"
else
    mongod
fi
