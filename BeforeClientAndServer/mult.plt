set terminal gif animate delay 100
set output 'MULT.gif'
stats 'mult.txt' nooutput

rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)

do for[i=1:int(STATS_blocks)]{
plot 'mult.txt' index(i-1) using 1:2:(rgb($3,$4,$5)) with points lc rgb variable
}
