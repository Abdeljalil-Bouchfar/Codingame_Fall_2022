#pragma GCC optimize("Ofast","unroll-loops", "omit-frame-pointer", "inline")
#pragma GCC option("arch=native", "tune=native", "no-zero-upper")
#pragma GCC target("rdrnd", "popcnt", "avx", "bmi2")

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define abs(x)  (x < 0) ? -x : x
#define index(x, y) (((y) * 24) + (x))


#define POPULATION_SIZE 30
#define SELECTION 10
#define ELITE 1
#define MAX_TIME 49
#define MAX_GENERATIONS 500000
#define CROSSOVER_RATE 1
#define MUTATION_RATE 0.5
#define BUILD_RATE 0.5
#define MATTER_USE_RATE 1
#define WIDTH 24
#define HEIGHT 12
#define HEIGHT_WIDTH 288
#define FARM_SCORE 20
#define FARM_MAX_LOSE 2
#define FIGHT_DES 2

int width, height;

struct Point {
    int x, y;
};

struct voronoi {
    int me, opp;
	bool valid;
};

struct Seed {
  Point point;
  char id;
};

struct Unit {
    Point pos;
    Point aiming;
	int des_to_edge;
};

struct Build_info {
    int pl = 1, ol = 0;
    int matter = 0;
	int opp_nighb_count = 0;
};

struct tile_data
{
    int scrap_amount;
    int owner;
    int units;
    int my_unit_count;
    int op_unit_count;
    int recycler;
    int can_build;
    int can_spawn;
    int in_range_of_recycler;
	bool is_edge_point = false;
};

struct player_data
{
	voronoi voronoi_state;
	int total_farm_recylers = 0;
	int matter;
	int initial_point;
	int recyclers_count;
	int cells_count;
	int edge_points_count;
	int build_cells_count;
	int spawn_cells_count;
	int units_count;
	int fight_units_count;
	int expand_units_count;
	Point *edge_points;
	Point *cells;
	Point *build_cells;
	Point *spawn_cells;
	Unit *units;
	Unit *fight_units;
	Unit *expand_units;
	// players map data: Positive number for player units and negative numbers for opponent units
	int *players_map_state;
	// Data of each point in the map from the player prospective
	tile_data *tiles;
};

struct game_data {
    int turn  = 0,phase = 0;
	player_data opp;
	player_data me;
	Seed *cell_seeds;
	int cell_seeds_count;
    // Here we have two maps [o_map] for real turn map
    char *tmap;
	// Voronoi diagram map
    char *vmap;
	// Voronoi diagram edge points
	// My units aiming
    map<pair<int, int>, vector<Point>> aiming;
    // This is a counter to track number of visited edge points bettwen turns
    int visited_edge_point_count = 0;
};

// Structure to represent an gene
struct Gene {
	char c;
	Point src;
	Point des;
};

// Structure to represent a chromosome
struct Chromosome {
	Gene mv_genes[100];
	Gene sb_genes[100];
	int mv_nbr;
	int sb_nbr;
	int score;
	int matter_cost;
	int info[10];
};

//---------------------- Data ------------------------

void set_plyaer_cells(game_data *gd, player_data *player)
{
	const int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	for(int i = 0; i < player->cells_count; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Point cell = player->cells[i];
			int x = cell.x + dx[j][0], y = cell.y + dx[j][1];
			if (x >= 0 && y >= 0 && x < width && y < height && player->tiles[index(x, y)].owner == 0 && gd->tmap[index(x, y)] != '#')
			{
				if (player->tiles[index(cell.x, cell.y)].units == 0 && player->tiles[index(cell.x, cell.y)].recycler == 0)
                	player->build_cells[player->build_cells_count++] = cell;
				if (player->tiles[index(cell.x, cell.y)].recycler == 0)
					player->spawn_cells[player->spawn_cells_count++] = cell;
				break;
			}
		}
	}
}

voronoi voronoi_state(game_data *gd)
{
    voronoi res;
    res.me = 0, res.opp = 0;
    char visited[HEIGHT_WIDTH] = {0};
    int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    int s = 0, e = gd->cell_seeds_count;
    Seed tmp;
    while (e > s)
    {
        for (int i = 0; i < 4; i++)
        {
            int x = gd->cell_seeds[s].point.x + dx[i][0], y = gd->cell_seeds[s].point.y + dx[i][1];
            if (x >= 0 && x < width && y >= 0 && y < height && !visited[index(x, y)] && gd->tmap[index(x, y)] != '#')
            {
                tmp.point.x = x, tmp.point.y = y, tmp.id = gd->cell_seeds[s].id;
                gd->cell_seeds[e++] = tmp, visited[index(x, y)] = gd->cell_seeds[s].id;
                res.me += tmp.id == 1, res.opp += tmp.id == 2;
            }
        }
        s++;
    }
    return res;
}


void set_data(game_data *gd)
{
	Point pos;
	Point t_aiming;
    scanf("%d%d", &(gd->me.matter), &(gd->opp.matter));
	gd->cell_seeds_count = 0;
    gd->me.recyclers_count = 0, gd->me.units_count = 0, gd->me.cells_count = 0, gd->me.build_cells_count = 0, gd->me.spawn_cells_count = 0;
    gd->opp.recyclers_count = 0, gd->opp.units_count = 0, gd->opp.cells_count = 0, gd->opp.build_cells_count = 0, gd->opp.spawn_cells_count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
			gd->me.tiles[index(j, i)].is_edge_point = false;
			if (gd->phase == 1)
			{
				gd->me.players_map_state[index(j, i)] = 0, gd->opp.players_map_state[index(j, i)] = 0;
				gd->vmap[index(j, i)] = 0;
			}
			pos.x = j, pos.y = i;
            scanf("%d%d%d%d%d%d%d", &(gd->me.tiles[index(j, i)].scrap_amount), &(gd->me.tiles[index(j, i)].owner), &(gd->me.tiles[index(j, i)].units), &(gd->me.tiles[index(j, i)].recycler), &(gd->me.tiles[index(j, i)].can_build), &(gd->me.tiles[index(j, i)].can_spawn), &(gd->me.tiles[index(j, i)].in_range_of_recycler));
			if (gd->me.tiles[index(j, i)].units)
			{
				gd->me.players_map_state[index(j, i)] = gd->me.tiles[index(j, i)].owner == 1 ? gd->me.tiles[index(j, i)].units : gd->me.tiles[index(j, i)].units * -1;
				gd->opp.players_map_state[index(j, i)] = gd->me.tiles[index(j, i)].owner == 0 ? gd->me.tiles[index(j, i)].units : gd->me.tiles[index(j, i)].units * -1;
			}


			// Save my units and openent units in two seperated vectors
            for (int u = 0; u < gd->me.tiles[index(j, i)].units; u++)
            {
				t_aiming = {-1, -1};
                if (gd->me.tiles[index(j, i)].owner == 1)
                {
					if (!gd->aiming[make_pair(j, i)].empty())
                    {
						// cerr << j << ", " << i << " -> " << gd->aiming[make_pair(j, i)][0].x << ", " << gd->aiming[make_pair(j, i)][0].y << endl;
                        t_aiming = *gd->aiming[make_pair(j, i)].begin();
                        gd->aiming[make_pair(j, i)].erase(gd->aiming[make_pair(j, i)].begin());
                    }
                    gd->me.units[gd->me.units_count++] = {pos, t_aiming, 0}, gd->me.tiles[index(j, i)].my_unit_count++;
                    gd->tmap[index(j, i)] = 'M';
                    // Here i set the my edge robots from up [mup] and down [mdp] 
                    if (gd->turn == 0 && i > 0 && gd->me.tiles[index(j, i - 1)].owner == 1)
                        gd->me.initial_point = j;
                }
                else
                {
                    gd->opp.units[gd->opp.units_count++] = {pos, t_aiming, 0}, gd->me.tiles[index(j, i)].op_unit_count++;
                    // add openent unit as seeds point for voronoi
                    gd->tmap[index(j, i)] = 'O';
					if (gd->turn == 0 && i > 0 && gd->me.tiles[index(j, i - 1)].owner == 0)
                        gd->opp.initial_point = j;
                }
            }
            // Count recyclers
            if (gd->me.tiles[index(j, i)].recycler && gd->phase == 0)
			{
                gd->me.total_farm_recylers += gd->me.tiles[index(j, i)].owner == 1;
				gd->opp.total_farm_recylers += gd->me.tiles[index(j, i)].owner == 0;
			}
            // Count cells and save my and opponenet cells points
            if (gd->me.tiles[index(j, i)].owner == 1)
			{
                gd->me.cells[gd->me.cells_count++] = pos;
				gd->cell_seeds[gd->cell_seeds_count].point = pos, gd->cell_seeds[gd->cell_seeds_count++].id = 1;

			}
            else if (gd->me.tiles[index(j, i)].owner == 0)
			{
                gd->opp.cells[gd->opp.cells_count++] = pos;
				gd->cell_seeds[gd->cell_seeds_count].point = pos, gd->cell_seeds[gd->cell_seeds_count++].id = 2;
			}

            // 2D map
            if (gd->me.tiles[index(j, i)].units == 0)
                gd->tmap[index(j, i)] = (gd->me.tiles[index(j, i)].scrap_amount == 0  || gd->me.tiles[index(j, i)].recycler) ? '#' : '.';
        }
    }
	voronoi vs = voronoi_state(gd);
	gd->me.voronoi_state = vs;
	int tmp = vs.me;
	vs.me = vs.opp, vs.opp = tmp;
	gd->opp.voronoi_state = vs;
	for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
    		int owner = gd->me.tiles[index(j, i)].owner;
			if (owner == 1)
				owner = 0;
			else if (owner == 0)
				owner = 1;
			gd->opp.tiles[index(j, i)].owner = owner;
			gd->opp.tiles[index(j, i)].units = gd->me.tiles[index(j, i)].units;
			gd->opp.tiles[index(j, i)].my_unit_count = gd->me.tiles[index(j, i)].op_unit_count ;
			gd->opp.tiles[index(j, i)].op_unit_count = gd->me.tiles[index(j, i)].my_unit_count ;
			gd->opp.tiles[index(j, i)].recycler = gd->me.tiles[index(j, i)].recycler ;
			gd->opp.tiles[index(j, i)].in_range_of_recycler = gd->me.tiles[index(j, i)].in_range_of_recycler ;
		}
	}
	set_plyaer_cells(gd, &gd->me);
	set_plyaer_cells(gd, &gd->opp);
}

//---------------------- Tools -----------------------

unsigned long	get_current_time(unsigned long st_time)
{
	struct timeval	tp;
	gettimeofday(&tp, NULL);
	return (((tp.tv_sec * 1000) + tp.tv_usec / 1000) - st_time);
}

int popu_comp(const void *p, const void *q) 
{
    return (((Chromosome *)q)->score - ((Chromosome *)p)->score);
}

void sort_units(Unit units[], int n) {
    for (int i = 1; i < n; i++) {
        Unit current = units[i];
        int j = i - 1;
        while (j >= 0 && units[j].des_to_edge > current.des_to_edge) {
            units[j + 1] = units[j];
            j--;
        }
        units[j + 1] = current;
    }
}

void print_popu_err(Chromosome *p, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p[i].mv_nbr; j++)
		{
			fprintf(stderr, "%c %i %i %i %i; ", p[i].mv_genes[j].c, p[i].mv_genes[j].src.x, p[i].mv_genes[j].src.y
			, p[i].mv_genes[j].des.x, p[i].mv_genes[j].des.y);
		}
		for (int j = 0; j < p[i].sb_nbr; j++)
			fprintf(stderr, "%c %i %i; ", p[i].sb_genes[j].c, p[i].sb_genes[j].des.x, p[i].sb_genes[j].des.y);
		fprintf(stderr, "| score: %i\n------------------------\n", p[i].score);
	}
}

void print_chrome(Chromosome chrome)
{

	for (int j = 0; j < chrome.mv_nbr; j++)
	{
		printf("MOVE 1 %i %i %i %i;", chrome.mv_genes[j].src.x, chrome.mv_genes[j].src.y
		, chrome.mv_genes[j].des.x, chrome.mv_genes[j].des.y);
	}
	for (int j = 0; j < chrome.sb_nbr; j++)
	{
		if (chrome.sb_genes[j].c == 'S')
			printf("SPAWN 1 %i %i;", chrome.sb_genes[j].des.x, chrome.sb_genes[j].des.y);
		else
			printf("BUILD %i %i;", chrome.sb_genes[j].des.x, chrome.sb_genes[j].des.y);
	}
}

void print_map(char *mess, char map[HEIGHT_WIDTH], game_data *gd)
{
	fprintf(stderr, "%s\n", mess);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
			fprintf(stderr, "%c ", map[index(x, y)]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "------------------------\n");
}

// ----------- Generate a random chromosome -------------------
Gene random_move(Point src, game_data *gd)
{
	Gene move;
	int c = 1;
	Point moves[5];
	moves[0] = src;
	Point p;
	int x = src.x, y = src.y;
	if (y > 0 && gd->tmap[index(x, y - 1)] != '#')
	{
		p.y = y - 1, p.x = x;
		moves[c++] = p;
	}
	if (x > 0 && gd->tmap[index(x - 1, y)] != '#')
	{
		p.y = y, p.x = x - 1;
		moves[c++] = p;
	}
	if (y + 1 < height && gd->tmap[index(x, y + 1)] != '#')
	{
		p.y = y + 1, p.x = x;
		moves[c++] = p;
	}
	if (x + 1 < width && gd->tmap[index(x + 1, y)] != '#')
	{
		p.y = y, p.x = x + 1;
		moves[c++] = p;
	}
	move.c = 'M';
	move.src = move.des = src;
	move.des = moves[rand() % c];
	return move;
}

void random_spawn_build(Gene *g, player_data *player, game_data *gd)
{
	if ((double)rand() / RAND_MAX < BUILD_RATE && player->build_cells_count)
	{
		g->c = 'B';
		g->src = g->des = player->build_cells[rand() % player->build_cells_count];
		gd->tmap[index(g->src.x, g->src.y)] = '#';
	}
	else if (player->spawn_cells_count)
	{
		g->src = g->des = player->spawn_cells[rand() % player->spawn_cells_count];
		for (int i = 0; gd->tmap[index(g->src.x, g->src.y)] == '#' && i < player->spawn_cells_count; i++)
			g->src = g->des = player->spawn_cells[rand() % player->spawn_cells_count];
		g->c = 'S';
		// gd->tmap[index(g->src.x, g->src.y)] = '#';
	}
}

void randomChromosome(Chromosome *s, player_data *player, game_data *gd)
{
	s->score = s->matter_cost = s->mv_nbr = s->sb_nbr = 0;
	// Pick a random spawn if matter and SPAWN_RATE
	for (int i = 0; player->matter - i >= 10; i += 10)
		if ((double)rand() / RAND_MAX < MATTER_USE_RATE)
		{
			Gene sbg;
			sbg.c = 0;
			random_spawn_build(&sbg, player, gd);
			if (sbg.c)
				s->sb_genes[s->sb_nbr++] = sbg;
		}
	
	// Pick a random move for each of unit
	for (int i = 0; i < player->fight_units_count; i++)
		s->mv_genes[s->mv_nbr++] = random_move(player->fight_units[i].pos, gd);

	for (int i = 0; i < s->sb_nbr; i++)
		if (s->sb_genes[i].c == 'B')
			gd->tmap[index(s->sb_genes[i].src.x, s->sb_genes[i].src.y)] = '.';
}

// Perform mutation on a chromosome at a random index
Chromosome mutate(Chromosome *s, game_data *gd, player_data *player)
{
	Chromosome mutated = *s;
	if (s->mv_nbr + s->sb_nbr)
	{
		int mutationIndex = rand() % (s->mv_nbr + s->sb_nbr);

		if (mutationIndex < s->mv_nbr)
			mutated.mv_genes[mutationIndex] = random_move(mutated.mv_genes[mutationIndex].src, gd);
		else
		{
			int i = s->mv_nbr + s->sb_nbr - mutationIndex - 1;
			random_spawn_build(&(mutated.sb_genes[i]), player, gd);
			if (mutated.sb_genes[i].c == 'B')
				gd->tmap[index(mutated.sb_genes[i].src.x, mutated.sb_genes[i].src.y)] = '.';
		}
	}
	return mutated;
}

// Perform crossover between two chromosomes at a random index
void crossover(Chromosome *child, Chromosome *s1, Chromosome *s2, game_data *gd)
{
	child->score = child->mv_nbr = child->matter_cost = child->sb_nbr = 0;
	// Moves crossover
	int crossoverIndex = 0;
	if (s1->mv_nbr && s2->mv_nbr)
	{
		crossoverIndex = rand() % min(s1->mv_nbr, s2->mv_nbr);
		for (int i = 0; i < crossoverIndex; i++)
			child->mv_genes[child->mv_nbr++] = s1->mv_genes[i];

		for (int i = crossoverIndex; i < s2->mv_nbr; i++)
			child->mv_genes[child->mv_nbr++] = s2->mv_genes[i];
	}

	// Spawn and build cross over
	if (s1->sb_nbr != 0 && s2->sb_nbr != 0)
		crossoverIndex = rand() % min(s1->sb_nbr, s2->sb_nbr);
	else
		crossoverIndex = s1->sb_nbr;
	for (int i = 0; i < crossoverIndex; i++)
		child->sb_genes[child->sb_nbr++] = s1->sb_genes[i];

	for (int i = crossoverIndex; i < s2->sb_nbr; i++)
		child->sb_genes[child->sb_nbr++] = s2->sb_genes[i];
}


void make(Chromosome *player, Chromosome *opp, game_data * gd, int *ge_map)
{
	// Sbawn and build
	// Mark Build points
	for (int i = 0; i < player->sb_nbr; i++)
	{
		int x = player->sb_genes[i].src.x, y = player->sb_genes[i].src.y;
		if (player->sb_genes[i].c == 'B')
		{
			ge_map[index(x, y)] = 1000;
			gd->tmap[index(x, y)] = '#';
		}
	}
	for (int i = 0; i < opp->sb_nbr; i++)
	{
		int x = opp->sb_genes[i].src.x, y = opp->sb_genes[i].src.y;
		if (opp->sb_genes[i].c == 'B')
		{
			ge_map[index(x, y)] = 1000;
			gd->tmap[index(x, y)] = '#';
		}
	}

	// make spawns
	for (int i = 0; i < player->sb_nbr; i++)
	{
		int x = player->sb_genes[i].src.x, y = player->sb_genes[i].src.y;
		if (player->sb_genes[i].c == 'S' && ge_map[index(x, y)] != 1000)
			ge_map[index(x, y)]++;
	}
	for (int i = 0; i < opp->sb_nbr; i++)
	{
		int x = opp->sb_genes[i].src.x, y = opp->sb_genes[i].src.y;
		if (opp->sb_genes[i].c == 'S' && ge_map[index(x, y)] != 1000)
			ge_map[index(x, y)]--;
	}

	// make Movment
	for (int i = 0; i < player->mv_nbr; i++)
	{
		if (ge_map[index(player->mv_genes[i].des.x, player->mv_genes[i].des.y)] != 1000)
		{
			int x = player->mv_genes[i].src.x, y = player->mv_genes[i].src.y;
			ge_map[index(x, y)]--;
			x = player->mv_genes[i].des.x, y = player->mv_genes[i].des.y;
			ge_map[index(x, y)]++;
		}
	}
	for (int i = 0; i < opp->mv_nbr; i++)
	{
		if (ge_map[index(opp->mv_genes[i].des.x, opp->mv_genes[i].des.y)] != 1000)
		{
			int x = opp->mv_genes[i].src.x, y = opp->mv_genes[i].src.y;
			ge_map[index(x, y)]++;
			x = opp->mv_genes[i].des.x, y = opp->mv_genes[i].des.y;
			ge_map[index(x, y)]--;
		}
	}
}

void undo(Chromosome *player, Chromosome *opp, game_data * gd, int *ge_map)
{
	// undo Movment
	for (int i = 0; i < player->mv_nbr; i++)
	{
		if (ge_map[index(player->mv_genes[i].des.x, player->mv_genes[i].des.y)] != 1000)
		{
			int x = player->mv_genes[i].src.x, y = player->mv_genes[i].src.y;
			ge_map[index(x, y)]++;
			x = player->mv_genes[i].des.x, y = player->mv_genes[i].des.y;
			ge_map[index(x, y)]--;
		}
	}
	for (int i = 0; i < opp->mv_nbr; i++)
	{
		if (ge_map[index(opp->mv_genes[i].des.x, opp->mv_genes[i].des.y)] != 1000)
		{
			int x = opp->mv_genes[i].src.x, y = opp->mv_genes[i].src.y;
			ge_map[index(x, y)]--;
			x = opp->mv_genes[i].des.x, y = opp->mv_genes[i].des.y;
			ge_map[index(x, y)]++;
		}
	}

	// undo SPAWN
	for (int i = 0; i < player->sb_nbr; i++)
	{
		int x = player->sb_genes[i].src.x, y = player->sb_genes[i].src.y;
		if (player->sb_genes[i].c == 'S' && ge_map[index(x, y)] != 1000)
			ge_map[index(x, y)]--;
	}
	for (int i = 0; i < opp->sb_nbr; i++)
	{
		int x = opp->sb_genes[i].src.x, y = opp->sb_genes[i].src.y;
		if (opp->sb_genes[i].c == 'S' && ge_map[index(x, y)] != 1000)
			ge_map[index(x, y)]++;
	}
	// undo Build
	for (int i = 0; i < player->sb_nbr; i++)
	{
		int x = player->sb_genes[i].src.x, y = player->sb_genes[i].src.y;
		if (ge_map[index(x, y)] == 1000)
		{
			ge_map[index(x, y)] = 0;
			gd->tmap[index(x, y)] = '.';
		}
	}
	for (int i = 0; i < opp->sb_nbr; i++)
	{
		int x = opp->sb_genes[i].src.x, y = opp->sb_genes[i].src.y;
		if (ge_map[index(x, y)] == 1000)
		{
			ge_map[index(x, y)] = 0;
			gd->tmap[index(x, y)] = '.';
		}
	}
} 

int gen_bfs(Point pos, game_data *gd, player_data *player, int depth)
{
	int des = 0;
	int qs = 0, qe = 0;
    Point que[300];
	que[qe++] = pos;
    int distances[HEIGHT_WIDTH] = {0};
    int visited[HEIGHT_WIDTH] = {0};
	visited[index(pos.x, pos.y)] = 1;
	int *state = player->players_map_state;
    while(qe > qs && des < depth)
    {
        int x = que[qs].x, y = que[qs++].y;
        des = distances[index(x, y)];
		if (state[index(x, y)] < 0 || player->tiles[index(x, y)].owner == 0)
			return des;
	
        if (x + 1 < width && gd->tmap[index(x + 1, y)] != '#' && !visited[index(x + 1, y)] && state[index(x + 1, y)] != 1000)
            que[qe++] = {x + 1, y}, visited[index(x + 1, y)] = 1, distances[index(x + 1, y)] += distances[index(x, y)] + 1;

        if (x > 0 && gd->tmap[index(x - 1, y)] != '#' && !visited[index(x - 1, y)] && state[index(x - 1, y)] != 1000)
            que[qe++] = {x - 1, y}, visited[index(x - 1, y)] = 1, distances[index(x - 1, y)] += distances[index(x, y)] + 1;
        
		if (y + 1 < height && gd->tmap[index(x, y + 1)] != '#' && !visited[index(x, y + 1)] && state[index(x, y + 1)] != 1000)
            que[qe++] = {x, y + 1}, visited[index(x, y + 1)] = 1, distances[index(x, y + 1)] += distances[index(x, y)] + 1;
        
        if (y > 0 && gd->tmap[index(x, y - 1)] != '#' && !visited[index(x, y - 1)] && state[index(x, y - 1)] != 1000)
            que[qe++] = {x, y - 1}, visited[index(x, y - 1)] = 1, distances[index(x, y - 1)] += distances[index(x, y)] + 1;

    }
    return 1000;
}

// Build scoring: oppenet unis count - (my losed tiles * 10) + (opponent losed tiles * 10) + total giving matter + matter state
Build_info get_build_info(Point pos, player_data *player, player_data *opp, game_data *gd)
{
	Build_info info;
	int tile_scm = player->tiles[index(pos.x, pos.y)].scrap_amount;
	const int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	for (int i = 0; i < 4; i++)
	{
		int x = pos.x + dx[i][0], y = pos.y + dx[i][1];
		if (x >= 0 && y >= 0 && x < width && y < height)
		{
			int tscrm = player->tiles[index(x, y)].scrap_amount;
			info.pl += (player->tiles[index(x, y)].owner == 1 && tscrm <= tile_scm);
			info.ol += (player->tiles[index(x, y)].owner == 0 && tscrm <= tile_scm);
			info.opp_nighb_count += player->tiles[index(x, y)].op_unit_count;
			if (tscrm >= tile_scm)
				info.matter += tile_scm;
			else
				info.matter += tscrm;
		}
	}
	return  info;
	//return 0;
}

// Calculate the fitness of a chromosome
// For scoring we gona have three vars for tiles progress: 
// - Controling Natural tiles -> %10
// - Controling Oppent tiles -> %40
// - Losing ower own tiles -> %50
// Spawn scoring: The given advantge of spawning by checking adjactent tiles state
int tg = 0;
void fitness(Chromosome *player_Chromosome, Chromosome *opp_Chromosome, game_data *gd, player_data *player, player_data *opp)
{
	player_Chromosome->score = 0;
	const int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	int nct = 0, oct = 0, plt = 0, o_nct = 0, adv_ps = 0, desadv_ps = 0;
	make(player_Chromosome, opp_Chromosome, gd, player->players_map_state);
	int owner, new_tile_val, tile_nighb_opp_count;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			owner = player->tiles[index(x, y)].owner;
			new_tile_val = player->players_map_state[index(x, y)];
			tile_nighb_opp_count = 0;
			for (int i = 0; i < 4; i++)
			{
				int tx = x + dx[i][0], ty = y + dx[i][1];
				if (x >= 0 && y >= 0 && x < width && y < height && player->players_map_state[index(tx, ty)] <= 0)
					tile_nighb_opp_count += abs(player->players_map_state[index(tx, ty)]);
			}
			adv_ps += (new_tile_val > tile_nighb_opp_count);
			desadv_ps += (new_tile_val < tile_nighb_opp_count);
			// Player Losed tile
			plt += (owner == 1 && new_tile_val < 0);
			// Opponent tile controlled
			oct += (owner == 0 && new_tile_val > 0);
			// Player controling natural tile
			nct += (owner == -1 && new_tile_val > 0);
			// Opponent controling natural tile
			o_nct += (owner == -1 && new_tile_val < 0);
		}
	}
	player_Chromosome->score += (plt * -100) + (nct * 50) + (oct * 60) + (o_nct * -60) + (adv_ps * 0) + (desadv_ps * 0);
	undo(player_Chromosome, opp_Chromosome, gd, player->players_map_state);
}

void genetic(game_data *gd)
{
	tg = 0;
	cerr << "genetic\n";
	const unsigned long	start_time = get_current_time(0);
	// Initialize the population with random chromosomes
	Chromosome *my_population = (Chromosome *) malloc(sizeof(Chromosome) * POPULATION_SIZE);
	Chromosome *op_population = (Chromosome *) malloc(sizeof(Chromosome) * POPULATION_SIZE);
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		randomChromosome(&(my_population[i]), &(gd->me), gd);
		randomChromosome(&(op_population[i]), &(gd->opp), gd);
	}
	int my_min, op_min;

	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		my_min = 100000, op_min = 100000;
		for (int j = 0; j < POPULATION_SIZE; j++)
		{
			fitness(&(my_population[i]), &(op_population[j]), gd, &gd->me, &gd->opp);
			if (my_population[i].score < my_min)
				my_min = my_population[i].score;
		
			fitness(&(op_population[i]), &(my_population[j]), gd, &gd->opp, &gd->me);
			if (op_population[i].score < op_min)
				op_min = op_population[i].score;
		}
		my_population[i].score = my_min;
		op_population[i].score = op_min;
	}
	qsort((void*)my_population, POPULATION_SIZE, sizeof(Chromosome), popu_comp);
	qsort((void*)op_population, POPULATION_SIZE, sizeof(Chromosome), popu_comp);
	

	// Run the genetic algorithm
	int g = 0;
	while (get_current_time(start_time) < MAX_TIME && g < MAX_GENERATIONS)
	{
		// Perform crossover and mutation
		for (int i = SELECTION; i < POPULATION_SIZE; i++)
		{
			crossover(&(my_population[i]), &(my_population[rand() % SELECTION]), &(my_population[rand() % SELECTION]), gd);
			crossover(&(op_population[i]), &(op_population[rand() % SELECTION]), &(op_population[rand() % SELECTION]), gd);
			if ((double)rand() / RAND_MAX < MUTATION_RATE)
			{
				my_population[i] = mutate(&(my_population[rand() % SELECTION]), gd, &(gd->me));
				op_population[i] = mutate(&(op_population[rand() % SELECTION]), gd, &(gd->opp));
			}
		}


		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			my_min = 100000, op_min = 100000;
			for (int j = 0; j < ELITE; j++)
			{
				fitness(&(my_population[i]), &(op_population[j]), gd, &gd->me, &gd->opp);
				if (my_population[i].score < my_min)
					my_min = my_population[i].score;
			
				fitness(&(op_population[i]), &(my_population[j]), gd, &gd->opp, &gd->me);
				if (op_population[i].score < op_min)
					op_min = op_population[i].score;
			}
			my_population[i].score = my_min;
			op_population[i].score = op_min;
		}
		qsort((void*)my_population, POPULATION_SIZE, sizeof(Chromosome), popu_comp);
		qsort((void*)op_population, POPULATION_SIZE, sizeof(Chromosome), popu_comp);
		g++;
		tg++;
		// print_popu_err(my_population, 1);
	}
	print_popu_err(my_population, 1);
	fprintf(stderr, "generations: %i\n", g);
	print_chrome(my_population[0]);
}
//----------------------- Expand (aka phase 0) --------------------------

int bfs_des = 0;
Point dynamic_bfs(Point pos, Point target, int target_key, game_data *gd, int depth)
{
    // target key is for what the bfs will search for
    // 0: specific target /\ 1: edge point in v_map /\ 2: for spawn cell,
	// 3: opp cells /\ 4: opp unit /\ 5: edge point in nodes
	// 6: opp cells || opp unit - 7: my cells || my units
	// if depth == 0 ignore it
    bfs_des = 0;
	int qs = 0, qe = 0;
    Point que[300];
	que[qe++] = pos;
    int distances[HEIGHT_WIDTH] = {0};
    int visited[HEIGHT_WIDTH] = {0};
	visited[index(pos.x, pos.y)] = 1;
    while(qe > qs)
    {
        int x = que[qs].x, y = que[qs++].y;
        bfs_des = distances[index(x, y)];

		if (depth && bfs_des == depth)	
			return {-1, -1};
        if (target_key == 0 && target.x == x && target.y == y)
            return {x, y};
        if (target_key == 1 && gd->vmap[index(x, y)] == '=')
            return {x, y};
        if (target_key == 2 && gd->me.tiles[index(x, y)].can_spawn)
            return {x, y};
        if (target_key == 3 && gd->me.tiles[index(x, y)].owner == 0)
            return {x, y};
        if (target_key == 4 && gd->tmap[index(x, y)] == 'O')
            return {x, y};
        if (target_key == 5 && gd->me.tiles[index(x, y)].is_edge_point)
            return {x, y};
        if (target_key == 6 && (gd->me.tiles[index(x, y)].op_unit_count || gd->me.tiles[index(x, y)].owner == 0))
            return {x, y};
        if (target_key == 7 && (gd->opp.tiles[index(x, y)].op_unit_count || gd->me.tiles[index(x, y)].owner == 1))
            return {x, y};


        if (x + 1 < width && gd->tmap[index(x + 1, y)] != '#' && !visited[index(x + 1, y)])
            que[qe++] = {x + 1, y}, visited[index(x + 1, y)] = 1, distances[index(x + 1, y)] += distances[index(x, y)] + 1;

        if (x > 0 && gd->tmap[index(x - 1, y)] != '#' && !visited[index(x - 1, y)])
            que[qe++] = {x - 1, y}, visited[index(x - 1, y)] = 1, distances[index(x - 1, y)] += distances[index(x, y)] + 1;
        
		if (y + 1 < height && gd->tmap[index(x, y + 1)] != '#' && !visited[index(x, y + 1)])
            que[qe++] = {x, y + 1}, visited[index(x, y + 1)] = 1, distances[index(x, y + 1)] += distances[index(x, y)] + 1;
        
        if (y > 0 && gd->tmap[index(x, y - 1)] != '#' && !visited[index(x, y - 1)])
            que[qe++] = {x, y - 1}, visited[index(x, y - 1)] = 1, distances[index(x, y - 1)] += distances[index(x, y)] + 1;

    }
    return {-1, -1};
}

// update this
void voronoi_line(game_data *gd)
{
	gd->me.edge_points_count = 0;
	const int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	int s = 0, e = 0;
	Seed seeds[600];
	Seed tmp;
	for (int i = 0; i < 4; i++)
	{
		tmp.id = 'M', tmp.point = gd->me.units[i].pos;
		seeds[e++] = tmp, gd->vmap[index(tmp.point.x, tmp.point.y)] = 'M';
		tmp.id = 'O', tmp.point = gd->opp.units[i].pos;
		seeds[e++] = tmp, gd->vmap[index(tmp.point.x, tmp.point.y)] = 'O';
	}

	
	while (e > s)
	{
		for (int i = 0; i < 4; i++)
		{
			int x = seeds[s].point.x + dx[i][0], y = seeds[s].point.y + dx[i][1];
			if (x >= 0 && x < width && y >= 0 && y < height && gd->tmap[index(x, y)] == '.')
			{
				if (gd->vmap[index(x, y)] && gd->vmap[index(x, y)] != seeds[s].id && gd->vmap[index(x, y)] == 'O')
				{
					gd->me.tiles[index(x, y)].is_edge_point = true;
					gd->vmap[index(x, y)] = '=';
					gd->me.edge_points[gd->me.edge_points_count++] = {x, y};
				}
				else if (!gd->vmap[index(x, y)])
				{
					tmp.point.x = x, tmp.point.y = y, tmp.id = seeds[s].id;
					seeds[e++] = tmp, gd->vmap[index(x, y)] = seeds[s].id;
				}
			}
		}
		s++;
	}
}

// Move to voronoi edge points
void voronoi_move(game_data *gd)
{
	for (int i = 0; i < gd->me.expand_units_count; i++)
	{
		dynamic_bfs(gd->me.expand_units[i].pos, {-1, -1}, 1, gd, 0);
		gd->me.expand_units[i].des_to_edge = bfs_des;
	}
	sort_units(gd->me.expand_units, gd->me.expand_units_count);

    for (int i = 0; i < gd->me.expand_units_count; i++)
    {
        if (gd->visited_edge_point_count >= gd->me.edge_points_count)
        {

            for (int i = 0; i < gd->me.edge_points_count; i++)
            	gd->vmap[index(gd->me.edge_points[i].x, gd->me.edge_points[i].y)] = '=';
            gd->visited_edge_point_count = 0;
        }
        Point pos = gd->me.expand_units[i].pos;
        Point edge_point;
		//
		if (gd->me.expand_units[i].aiming.x == -1)
			edge_point = dynamic_bfs(gd->me.expand_units[i].pos, gd->me.expand_units[i].aiming, 1, gd, 0);
		else
			edge_point = dynamic_bfs(gd->me.expand_units[i].pos, gd->me.expand_units[i].aiming, 0, gd, 0);
		//
		Point move = edge_point;
        int left_des = 100000, right_des = 100000, up_des = 100000, down_des = 100000;
        if (edge_point.x != -1)
        {
            if (pos.y + 1 < height && gd->tmap[index(pos.x, pos.y + 1)] != '#')
            {
                Point t = dynamic_bfs({pos.x, pos.y + 1}, edge_point, 0, gd, 0);
                if (t.y != -1)
                    down_des = bfs_des;
            }
            if (pos.y > 0 && gd->tmap[index(pos.x, pos.y - 1)] != '#')
            {
                Point t = dynamic_bfs({pos.x, pos.y - 1}, edge_point, 0, gd, 0);
                if (t.y != -1)
                    up_des = bfs_des;
            }
            if (pos.x + 1 < width && gd->tmap[index(pos.x + 1, pos.y)] != '#')
            {
                Point t = dynamic_bfs({pos.x + 1, pos.y}, edge_point, 0, gd, 0);
                if (t.y != -1)
                    right_des = bfs_des;
            }
            if (pos.x > 0 && gd->tmap[index(pos.x - 1, pos.y)] != '#')
            {
                Point t = dynamic_bfs({pos.x - 1, pos.y}, edge_point, 0, gd, 0);
                if (t.y != -1)
                    left_des = bfs_des;
            }
            if (min(down_des, min(up_des, min(right_des, left_des))) == down_des)
                move = {pos.x, pos.y + 1};
            else if (min(down_des, min(up_des, min(right_des, left_des))) == up_des)
                move = {pos.x, pos.y - 1};
            else if (min(down_des, min(up_des, min(right_des, left_des))) == right_des)
                move = {pos.x + 1, pos.y};
            else if (min(down_des, min(up_des, min(right_des, left_des))) == left_des)
                move = {pos.x - 1, pos.y};
            cout << "MOVE 1 " << pos.x << ' ' << pos.y << ' ' << move.x << ' ' << move.y << ';';
            gd->aiming[make_pair(move.x, move.y)].push_back(edge_point);
            gd->visited_edge_point_count += gd->vmap[index(edge_point.x, edge_point.y)] == '=';
            gd->vmap[index(edge_point.x, edge_point.y)] = '.';
			gd->me.tiles[index(move.x, move.y)].owner = 1;
			gd->me.tiles[index(move.x, move.y)].units++;
			gd->me.tiles[index(move.x, move.y)].my_unit_count++;
			gd->me.tiles[index(move.x, move.y)].can_build = 0;
			gd->tmap[index(move.x, move.y)] = 'M';
			gd->me.players_map_state[index(pos.x, pos.y)]--;
			gd->me.players_map_state[index(move.x, move.y)]++;
        }
        else
            cout << "WAIT;";
    }
}

void expand_spawn(game_data *gd)
{
	for (int i = 0; i < gd->me.edge_points_count; i++)
    {
		if (gd->me.matter < 10)
			break;
        Point c = gd->me.edge_points[i];
        if (gd->vmap[index(c.x, c.y)] == '=')
        {
        	Point res = dynamic_bfs(c, {-1, -1}, 2, gd, 0);
            if (res.x != -1)
            {
                cout << "SPAWN 1 " << res.x << ' ' << res.y << ";";
				gd->me.tiles[index(res.x, res.y)].can_build = 0;
				gd->tmap[index(res.x, res.y)] = 'M';
				gd->me.players_map_state[index(res.x, res.y)]++;
                gd->me.matter -= 10;
            }
        }
    }
}


void collect(game_data *gd)
{
	cerr << "Collect\n";
	const int dx[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	for (int j = 0; j < gd->me.cells_count && gd->me.matter >= 10; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			int x = gd->me.cells[j].x + dx[i][0], y = gd->me.cells[j].y + dx[i][1];
			if (x >= 0 && y >= 0 && x < width && y < height && gd->tmap[index(x, y)] != '#' && gd->me.tiles[index(x, y)].owner == -1)
			{
				cout << "SPAWN 1 " << gd->me.cells[j].x << ' ' << gd->me.cells[j].y << ";";
				gd->me.matter -= 10;
				gd->me.tiles[index(gd->me.cells[j].x, gd->me.cells[j].y)].my_unit_count++;
				break;
			}
		}
	}
	for (int j = 0; j < gd->me.units_count; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			int x = gd->me.units[j].pos.x + dx[i][0], y = gd->me.units[j].pos.y + dx[i][1];
			if (x >= 0 && y >= 0 && x < width && y < height && gd->me.tiles[index(x, y)].owner == -1 && gd->tmap[index(x, y)] != '#')
			{
				cout << "MOVE 1 " << gd->me.units[j].pos.x << ' ' << gd->me.units[j].pos.y << ' ' << x << ' ' << y << ";";
				break;
			}
		}
	}
}

void print_tiles_state(game_data *gd)
{

	cerr << "    ";
	for (int i = 0; i < width; i++)
	{
		if (i < 10)
			cerr << " ";
		cerr << " " << i;
	}
	cerr << endl;
	for (int i = 0; i < height; i++)
	{
		cerr << " " << i << " | ";
		for (int j = 0; j < width; j++)
		{
			if (gd->me.tiles[index(j, i)].owner >= 0)
				cerr << " " << gd->me.tiles[index(j, i)].owner << " ";
			else
				cerr << gd->me.tiles[index(j, i)].owner << " ";
		}
		cerr << "\n";
	}
	cerr << "-----------------------\n";
}

void farm_build(game_data *gd)
{
	vector<Point> lcp;
    cerr << "farm_build\n";
    int max_score = FARM_SCORE, x = -1, y = -1, scrm = 0;
    for (int c = 0; c < gd->me.cells_count; c++) {
		lcp = {};
        int tmp_score = 0, lose = 1, i = gd->me.cells[c].y, j = gd->me.cells[c].x;
		lcp.push_back({j , i});
        tmp_score = scrm = gd->me.tiles[index(j, i)].scrap_amount;
        if (gd->me.tiles[index(j, i)].can_build && !gd->me.tiles[index(j, i)].in_range_of_recycler)
        {
            if (j > 0 && gd->me.tiles[index(j - 1, i)].scrap_amount > 0)
            {
                if (gd->me.tiles[index(j - 1, i)].scrap_amount <= scrm)
                    tmp_score += gd->me.tiles[index(j - 1, i)].scrap_amount, lose += gd->me.tiles[index(j - 1, i)].owner == 1, lcp.push_back({j - 1, i});
                else
                    tmp_score += scrm;
            }
            
            if (i > 0 && gd->me.tiles[index(j ,i - 1)].scrap_amount > 0)
            {
                if (gd->me.tiles[index(j ,i - 1)].scrap_amount <= scrm)
                    tmp_score += gd->me.tiles[index(j ,i - 1)].scrap_amount, lose += gd->me.tiles[index(j ,i - 1)].owner == 1, lcp.push_back({j, i - 1});
                else
                    tmp_score += scrm;
            }
            
            if (j + 1 < width && gd->me.tiles[index(j + 1, i)].scrap_amount > 0)
            {
                if (gd->me.tiles[index(j + 1, i)].scrap_amount <= scrm)
                    tmp_score += gd->me.tiles[index(j + 1, i)].scrap_amount, lose += gd->me.tiles[index(j + 1, i)].owner == 1, lcp.push_back({j + 1, i});
                else
                    tmp_score += scrm;
            }
            
            if (i + 1 < height && gd->me.tiles[index(j, i + 1)].scrap_amount > 0)
            {
                if (gd->me.tiles[index(j, i + 1)].scrap_amount <= scrm)
                    tmp_score += gd->me.tiles[index(j, i + 1)].scrap_amount, lose += gd->me.tiles[index(j, i + 1)].owner == 1, lcp.push_back({j, i + 1});
                else
                    tmp_score += scrm;
            }
			for (auto &p : lcp)
				gd->tmap[index(p.x, p.y)] = '#';
			voronoi ts = voronoi_state(gd);
            if (tmp_score >= max_score && lose <= FARM_MAX_LOSE && gd->me.voronoi_state.me - ts.me <= lose)
                max_score = tmp_score, x = j, y = i;
			for (auto &p : lcp)
				gd->tmap[index(p.x, p.y)] = '.';            
        }
    }
    if (x != -1)
    {  
        cout << "BUILD " << x << ' ' << y << ";";
        gd->me.tiles[index(x, y)].recycler = 1;
        gd->me.matter -= 10;
    }
}

void update_phase(game_data *gd)
{
	if (gd->phase == 0)
	{
		for (int i = 0; i < gd->me.units_count; i++)
		{
			Point res = dynamic_bfs(gd->me.units[i].pos, {-1, -1}, 6, gd, 2);
			if (res.x != -1)
			{
				gd->phase = 1;
				return;
			}
		}
	}
	else if (gd->phase == 1)
	{
		for (int i = 0; i < gd->me.cells_count; i++)
    	{
			Point c = gd->me.cells[i];
    	    if (gd->tmap[index(c.x, c.y)] != '#')
    	    {
    	        Point res = dynamic_bfs(c, {-1, -1}, 3, gd, 0);
    	        if (res.x != -1)
    	            return;
    	    }
    	}
		gd->phase = 2;
	}
}

void routing_units(game_data *gd)
{
	gd->me.fight_units_count = gd->me.expand_units_count = 0;
	for (int i = 0; i < gd->me.units_count; i++)
	{
		Point res = dynamic_bfs(gd->me.units[i].pos, {-1, -1}, 6, gd, 3);
		if (res.x != -1)
			gd->me.fight_units[gd->me.fight_units_count++] = gd->me.units[i];
		else
			gd->me.expand_units[gd->me.expand_units_count++] = gd->me.units[i];
	}
	gd->opp.fight_units_count = gd->opp.expand_units_count = 0;
	for (int i = 0; i < gd->opp.units_count; i++)
	{
		Point res = dynamic_bfs(gd->opp.units[i].pos, {-1, -1}, 7, gd, 3);
		if (res.x != -1)
			gd->opp.fight_units[gd->opp.fight_units_count++] = gd->opp.units[i];
		else
			gd->opp.expand_units[gd->opp.expand_units_count++] = gd->opp.units[i];
	}
}

int main()
{
	freopen("/dev/null", "w", stderr);
	game_data gd;
	srand(time(NULL));  // Initialize the random number generator with the current time
    scanf("%d%d", &width, &height);
	gd.turn = gd.me.edge_points_count = gd.opp.edge_points_count = 0;
	gd.me.total_farm_recylers = gd.opp.total_farm_recylers = 0;
	gd.me.units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.me.fight_units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.me.expand_units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.me.cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.me.spawn_cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.me.build_cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.me.edge_points = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.me.players_map_state = (int *) malloc(sizeof(int) * HEIGHT_WIDTH);
	gd.me.tiles = (tile_data *) malloc(sizeof(tile_data) * HEIGHT_WIDTH);
	//
	gd.opp.units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.opp.fight_units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.opp.expand_units = (Unit *) malloc(sizeof(Unit) * HEIGHT_WIDTH);
	gd.opp.cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.opp.spawn_cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.opp.build_cells = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.opp.edge_points = (Point *) malloc(sizeof(Point) * HEIGHT_WIDTH);
	gd.opp.players_map_state = (int *) malloc(sizeof(int) * HEIGHT_WIDTH);
	gd.opp.tiles = (tile_data *) malloc(sizeof(tile_data) * HEIGHT_WIDTH);
	//
	gd.vmap = (char *) calloc(sizeof(char), HEIGHT_WIDTH);
	gd.tmap = (char *) calloc(sizeof(char), HEIGHT_WIDTH);
	gd.cell_seeds = (Seed *) malloc(sizeof(Seed) * HEIGHT_WIDTH * 10);


	bool flg = false;

    // game loop
    while (1) {
		set_data(&gd);
		if (gd.turn == 0 || gd.phase == 1)
			voronoi_line(&gd);
		update_phase(&gd);
		print_map("voronoi map: ", gd.vmap, &gd);
		if (gd.phase != 2)
		{
			routing_units(&gd);
			if (gd.me.matter >= 10 && ((gd.phase == 0 && gd.me.total_farm_recylers <= (width * height) / 30) || gd.me.total_farm_recylers <= gd.opp.total_farm_recylers))
				farm_build(&gd);
			if (gd.phase == 0)
				expand_spawn(&gd);
			voronoi_move(&gd);
			genetic(&gd);
		} 
		else
			collect(&gd);
        // To debug: fprintf(stderr, "Debug messages...\n");
		// print_tiles_state(&gd);
		printf("WAIT;");
		printf("MESSAGE P%i | M%i | O%i;\n", gd.phase, gd.me.voronoi_state.me, gd.me.voronoi_state.opp);
		gd.turn++;
    }
    return 0;
}
