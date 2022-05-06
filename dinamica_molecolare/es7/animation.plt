set datafile sep '\t'
set grid

set term qt pos 0,300
set view equal xyz
padding = L/2
set xrange [-padding:L+padding]
set yrange [-padding:L+padding]
set zrange [-padding:L+padding]

set style fill transparent solid 0.5 border border lc rgb 'dark-green'
set object 1 polygon from 0,0,0 to L,0,0 to L,L,0 to 0,L,0 to 0,0,0 fillcolor rgb 'light-green'
set object 2 polygon from 0,0,0 to L,0,0 to L,0,L to 0,0,L to 0,0,0 fillcolor rgb 'light-green'
set object 3 polygon from 0,0,0 to 0,L,0 to 0,L,L to 0,0,L to 0,0,0 fillcolor rgb 'light-green'
set object 4 polygon from L,0,0 to L,L,0 to L,L,L to L,0,L to L,0,0 fillcolor rgb 'light-green'
set object 5 polygon from L,L,0 to 0,L,0 to 0,L,L to L,L,L to L,L,0 fillcolor rgb 'light-green'
set object 6 polygon from 0,0,L to L,0,L to L,L,L to 0,L,L to 0,0,L fillcolor rgb 'light-green'

n = 0
do for [step=1:N_step] {
	#set output sprintf('png/animation%03.0f.png',step)
	splot "coordinate.xyz" every ::n::n+N u 2:3:4 pt 7 ps 0.5 lc rgb "red" title sprintf('step = %03.0f    t = %f',step,n*dt/N)#,\
		  #"coordinate.xyz" every ::n::n+N u 2:3:4:5:6:7 with vectors filled head lw 2 lc rgb "blue" title ""
	      #"coordinate.xyz" every ::n::n+N u 2:3:4:1 w labels offset 2 title ""
	
	n=n+N*skip
	pause pausa
}


pause-1