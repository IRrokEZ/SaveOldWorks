rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) 
plot 'MyEmTest1.txt' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable