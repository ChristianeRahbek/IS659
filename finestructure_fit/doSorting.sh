for run in 165 166; do
    SETUP="../setup/setup.json"

    echo "Sorting run : $run using $SETUP"
    FILES="../runs/Run${run}*"
    echo "$FILES"

    Sorter --ignore C1,C2,C3,C4,P1,P2,P3,P4 -s $SETUP -m ../sorting/matcher.json $FILES
done

