# Game Fall Challenge 2022

## STATEMENT
Each player controls a team of robots.

Robots are deployed in a field of abandoned electronics, their purpose is to refurbish patches of this field into functional tech.

The robots are also capable of self-disassembly and self-replication, but they need raw materials from structures called Recyclers which the robots can build.

The structures will recycle everything around them into raw matter, essentially removing the patches of electronics and revealing the Grass below.

Players control a team of these robots in the midst of a playful competition to see which team can control the most patches of a given scrap field. They do so by marking patches with their team's color, all with the following constraints:

If robots of both teams end up on the same patch, they must disassemble themselves one for one. The robots are therefore removed from the game, only leaving at most one team on that patch.
The robots may not cross the grass, robots that are still on a patch when it is completely recycled must therefore disassemble themselves too.
Once the games are over, the robots will dutifully re-assemble and go back to work as normal.

## The game looks like this:
<img src="https://github.com/Abdeljalil-Bouchfar/Codingame_Fall_2022/blob/master/imgs/example.png?raw=true" alt="">

## My approach
In the beginning and because of the perfect input information I thought that going with a heuristic approach would work for this game, but when leagues progressed and the code got complicated more and more I found myself unable to make any improvements.

Then one of my colleagues (Sheeesh--- 👨🏻‍💻) suggested using a Genetic algorithm with two populations of individuals (using 30 in size worked perfectly for me) one for me and one for the opponent, and the best 10 individuals will be progressed to the next generation, and the next generation will be played against those best individuals. That way, both my population and the opponent's population will keep improving until I run out of time (50ms)
Finally, I play the best individual in my population.

### My overall ranking
<img src="https://github.com/Abdeljalil-Bouchfar/Codingame_Fall_2022/blob/master/imgs/indiv_rand.png?raw=true" alt="">

### My school rank (1337 ❤️)
<img src="https://github.com/Abdeljalil-Bouchfar/Codingame_Fall_2022/blob/master/imgs/team_rand.png?raw=true" alt="">
