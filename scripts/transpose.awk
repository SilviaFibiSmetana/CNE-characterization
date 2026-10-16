BEGIN { FS = OFS = "\t" }
{
    # For each line, store each field in the respective array element
    for (i = 1; i <= NF; i++) {
        row[i] = (row[i] ? row[i] OFS : "") $i
    }
}
ENDFILE {
    # Print each row and reset for the next file
    for (i = 1; row[i] != ""; i++) {
        print row[i]
        delete row[i]
    }
}

