for run in 121; do
    SETUP="../setup/setup.json"

    echo "Sorting run : $run using $SETUP"
    FILES="../runs/Run${run}*"
    echo "$FILES"

    Sorter --ignore C1,C2,C3,C4 -s $SETUP -m ../sorting/matcher.json $FILES
done

