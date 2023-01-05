#include "finders.h"
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdbool.h>

#define MC MC_1_17
// how many threads to use (unless overwritten)
#define THREADS 4
// how many strongholds to check for each seed
#define STRONGHOLDS 50
// these are in chunks
#define MUSH_TEST 7
#define SOUL_TEST 7
// these are in blocks, but will be converted to chunks and/or regions
#define FORTRESS_DIST 1000
#define OUTPOST_DIST 10000
#define HUTS_DIST 20000
#define MONUMENTS_DIST 20000
#define BIOMES_DIST 20000
#define JUNGLE_TEST 1000
#define TAIGA_TEST 1000
#define QUARRY_LENGTH 1000
#define SLIMES_DIST 80000
// how many slime chunks you want within the perimeter
#define SLIMES_MIN 40
// making this higher will be very slightly faster, but might miss some slime locations
#define SLIMES_CUTOFF (SLIMES_MIN / 4)

static int distpos(Pos, Pos);

static int count_hut;
static int count_monum;
static StructureConfig conf_fortress;
static StructureConfig conf_outpost;
static StructureConfig conf_hut;
static StructureConfig conf_monument;

static inline int *reallocCache(int *cache, size_t *current, const Generator *g, Range r)
{
	size_t len = getMinCacheSize(g, r.scale, r.sx, r.sy, r.sz);
	
	if(len > *current)
	{
		cache = (int *)realloc(cache, len * sizeof(int));
		*current = len;
	}
	
	return cache;
}

// this function is a lot of copy-paste, so gonna make some temp. macros
static inline bool haveMaxEnclosure(Pos p1, Pos p2, Pos p3, Pos p4, int radius)
{
	#define testdist(a,b) if(distpos(a, b) > 2*radius) { return false; }
	testdist(p1, p2)
	testdist(p1, p3)
	testdist(p2, p3)
	testdist(p1, p4)
	testdist(p2, p4)
	testdist(p3, p4)
	#undef testdist
	
    int x1 = p1.x, z1 = p1.z;
	int x2 = p2.x, z2 = p2.z;
	int x3 = p3.x, z3 = p3.z;
	int x4 = p4.x, z4 = p4.z;

    // sort x and z values
	#define sort(a,b) if(a > b) { int h = a; a = b; b = h; }
	sort(x1, x2) sort(x2, x3) sort(x3, x4)
	sort(x1, x2) sort(x2, x3)
	sort(x1, x2)
	sort(z1, z2) sort(z2, z3) sort(z3, z4)
	sort(z1, z2) sort(z2, z3)
	sort(z1, z2)
	#undef sort
	
    int minx = x4 - radius;
    int maxx = x1 + radius;
    int minz = z4 - radius;
    int maxz = z1 + radius;

    // brute force the best center position
    for(int z = minz; z <= maxz; z++)
        for(int x = minx; x <= maxx; x++)
        {
            int sq = (x-x1)*(x-x1) + (z-z1)*(z-z1);
			sq = max(sq, (x-x2)*(x-x2) + (z-z2)*(z-z2));
			sq = max(sq, (x-x3)*(x-x3) + (z-z3)*(z-z3));
			sq = max(sq, (x-x4)*(x-x4) + (z-z4)*(z-z4));
            if(sq <= radius * radius)
				return true;
        }

    return false;
}

static bool isMushroom(Generator *g, Pos pos)
{
	static __thread int *biomeIds = NULL;
	static __thread size_t idsStorage = 0;
	
	Range r = {16, pos.x / 16 - MUSH_TEST / 2, pos.z / 16 - MUSH_TEST / 2, MUSH_TEST, MUSH_TEST, 15, 1};
	biomeIds = reallocCache(biomeIds, &idsStorage, g, r);
	
	for(int dist = 1; dist <= MUSH_TEST; dist += 2)
	{
		Range r = {16, pos.x / 16 - dist / 2, pos.z / 16 - dist / 2, dist, dist, 15, 1};

		genBiomes(g, biomeIds, r);

		bool valid = true;
		for(int i = 0; i < dist * dist; ++i)
			valid &= biomeIds[i] == mushroom_fields;
		if(!valid)
			return false;
	}

	return true;
}

static bool hasFortress(Generator *g, int x, int z, Pos *ftpos)
{
	int span = 16 * conf_fortress.regionSize;
	
	int fdist = FORTRESS_DIST / span;
	int rx = x / span;
	int rz = z / span;

	Range range = {16, 0, 0, SOUL_TEST, SOUL_TEST, 15, 1};
	
	int posx = 0, posz = 0, dirx = 0, dirz = -1;
	for(int c = 0; c < (2 * fdist + 1) * (2 * fdist + 1); ++c)
	{
		int fortress = getStructurePos(Fortress, g->mc, g->seed, rx + posx, rz + posz, ftpos);
		if(fortress) {
			range.x = ftpos->x / 16 - SOUL_TEST / 2;
			range.z = ftpos->z / 16 - SOUL_TEST / 2;

			int *biomeIds = allocCache(g, range);
			genBiomes(g, biomeIds, range);

			bool valid = true;
			for(int i = 0; i < SOUL_TEST * SOUL_TEST; ++i)
				valid &= biomeIds[i] == soul_sand_valley;
			
			free(biomeIds);
			if(valid) return true;
		}
		
		posx += dirx;
		posz += dirz;
		if(dirz == -1 ? (posx == posz + 1) : (abs(posx) == abs(posz)))
		{
			int ndirx = -dirz;
			dirz = dirx;
			dirx = ndirx;
		}
	}

	return false;
}

static bool hasOutpost(Generator *g, Pos pos, Pos *oppos)
{
	int span = 16 * conf_outpost.regionSize;

	int rdist = OUTPOST_DIST / span;
	int rx = pos.x / span;
	int rz = pos.z / span;
	
	int x = 0, z = 0, dirx = 0, dirz = -1;
	for(int c = 0; c < (2 * rdist + 1) * (2 * rdist + 1); ++c)
	{
		int outpost = getStructurePos(Outpost, g->mc, g->seed, rx + x, rz + z, oppos);
		if(outpost)
			return true;
		
		x += dirx;
		z += dirz;
		if(dirz == -1 ? (x == z + 1) : (abs(x) == abs(z)))
		{
			int ndirx = -dirz;
			dirz = dirx;
			dirx = ndirx;
		}
	}

	return false;
}

static int isJungle(int biome)
{
	return (biome & ~0x80) == jungle || (biome & ~0x80) == jungle_hills || biome == jungle_edge || (biome & ~1) == bamboo_jungle;
}

static int isTaiga(int id)
{
	switch(id)
	{
	case giant_tree_taiga:
	case giant_tree_taiga_hills:
	case giant_spruce_taiga:
	case giant_spruce_taiga_hills:
		return 1;
	default:
		return 0;
	}
}

static bool hasBiome(int *biomeIds, int span, int testw, int testh, int (*biome)(int), Pos *pos)
{
	for(int z = 0; z + testh <= span; ++z)
		for(int x = 0; x + testw <= span; ++x)
		{
			int index = z * span + x;
			bool test = (*biome)(biomeIds[index])
					&& (*biome)(biomeIds[index + testw - 1])
					&& (*biome)(biomeIds[index + span * (testh - 1)])
					&& (*biome)(biomeIds[index + testw - 1 + span * (testh - 1)]);
			if(test)
			{
				bool valid = true;
				for(int dz = z; dz < z + testh; ++dz)
					for(int dx = x; dx < x + testw; ++dx)
						valid &= (*biome)(biomeIds[span * dz + dx]);

				if(valid)
				{
					pos->x = x;
					pos->z = z;
					return true;
				}
			}
		}
	return false;
}

// oppos: output positions
//  1. Jungle
//  2. Taiga
//  3. Mesa
//  4. Warm Ocean (not implemented)
//  5. Frozen Ocean (not implemented)
//  6. Desert
static bool hasBiomes(Generator *g, Pos pos, Pos *oppos)
{
	static __thread int *biomeIds = NULL;
	static __thread size_t idsStorage = 0;
	
	int span = BIOMES_DIST / 16; 

	Range r = {16, pos.x / 16 - span, pos.z / 16 - span, 2 * span + 1, 2 * span + 1, 15, 1};
	biomeIds = reallocCache(biomeIds, &idsStorage, g, r);
	genBiomes(g, biomeIds, r);

	bool valid = true
			&& hasBiome(biomeIds, 2 * span + 1, JUNGLE_TEST / 16, JUNGLE_TEST / 16, &isJungle, oppos)
			&& hasBiome(biomeIds, 2 * span + 1, TAIGA_TEST / 16, TAIGA_TEST / 16, &isTaiga, oppos + 1)
			&& (hasBiome(biomeIds, 2 * span + 1, QUARRY_LENGTH / 16, 2, &isMesa, oppos + 2)
				|| hasBiome(biomeIds, 2 * span + 1, 2, QUARRY_LENGTH / 16, &isMesa, oppos + 2))
		;
	
	for(int i = 0; i < 3; ++i)
	{
		oppos[i].x = 16 * oppos[i].x + pos.x - 8 * span;
		oppos[i].z = 16 * oppos[i].z + pos.z - 8 * span;
	}
	
	// use spiral code for finding _closest_ ocean/desert
	int x = 0, z = 0, dirx = 0, dirz = -1;
	bool desert_found = false;
	for(int c = 0; c < span * span && !desert_found; ++c)
	{
		int biome = biomeIds[(z + span/2) * span + (x + span/2 - 1)];
		if(!desert_found && biome == desert) {
			desert_found = true;
			oppos[5].x = pos.x + 16 * (x - 1);
			oppos[5].z = pos.z + 16 * z;
		}
		
		x += dirx;
		z += dirz;
		if(dirz == -1 ? (x == z + 1) : (abs(x) == abs(z)))
		{
			int ndirx = -dirz;
			dirz = dirx;
			dirx = ndirx;
		}
	}
	if(!desert_found) oppos[5].x = oppos[5].z = 100 * BIOMES_DIST - 1;

	return valid;
}

static int countMushs(Generator *g, Pos pos)
{
	Range r = {1, pos.x - 256, pos.z - 256, 512, 512, 63, 1};
	
	int *biomeIds = allocCache(g, r);
	genBiomes(g, biomeIds, r);
	
	int count = 0;
	for(int ix = 0; ix < 512 * 512; ++ix)
		count += (biomeIds[ix] == mushroom_fields);
	
	free(biomeIds);
	return count;
}

static const int secants[] = {0, 3, 5, 6, 6, 7, 7, 7, 8, 7, 7, 7, 6, 6, 5, 3, 0};
static int findSlimes(uint64_t seed, Pos starting_pos, Pos *slime_pos)
{
	bool slimes[24 * 32];
	
	// first round of testing is done in regions 8 x 8 chunks
	// these are chunk coordinates
	int baseX = starting_pos.x / 16 - 4, baseZ = starting_pos.z / 16 - 4;
	int lumpX = 0, lumpZ = 0;
	int dlx = 0, dlz = -8;
	
	int span = SLIMES_DIST / 16 / 8;
	for(int c = 0; c < (2 * span + 1) * (2 * span + 1); ++c)
	{
		// test lump
		int score = 0;
		for(int x = 0; x < 8; ++x)
		for(int z = 0; z < 8; ++z)
			score += isSlimeChunk(seed, baseX + lumpX + x, baseZ + lumpZ + z);
		
		// this rules out about 98% of non-viable chunks
		if(score >= SLIMES_CUTOFF)
		{
			// so we don't need to recalculate this 12608 times
			for(int dx = 0; dx < 24; ++dx)
			for(int dz = 0; dz < 24; ++dz)
				slimes[32 * dz + dx] = isSlimeChunk(seed, baseX + lumpX + dx - 8, baseZ + lumpZ + dz - 8);
			
			bool found = false;
			int centerX = 8;
			int centerZ = 8;
			int count;
			for(; !found && centerX < 16; ++centerX)
			for(; !found && centerZ < 16; ++centerZ)
			{
				count = 0;
				for(int dx = -8; dx <= 8; ++dx) for(int dz = -secants[dx+8]; dz <= secants[dx+8]; ++dz)
					count += slimes[32 * (centerZ + dz) + (centerX + dx)];
				if(count >= SLIMES_MIN)
					found = true;
			}
			
			if(found)
			{
				int chunkX = lumpX + centerX - 9;
				int chunkZ = lumpZ + centerZ - 9;
				slime_pos->x = 16 * chunkX + 8;
				slime_pos->z = 16 * chunkZ + 8;
				return count;
			}
		}
		
		// advance in spiral
		lumpX += dlx;
		lumpZ += dlz;
		if(dlz == -8 ? (lumpZ == lumpX-8) : (abs(lumpX) == abs(lumpZ)))
		{
			int ndlx = -dlz;
			int ndlz = dlx;
			dlx = ndlx;
			dlz = ndlz;
		}
	}
	
	return 0;
}

static bool hasRegDoubleStructure(enum StructureType type, Generator *g, int regX, int regZ, int radius, Pos *oppos)
{
	Pos mpos1, mpos2;
	
	int monum1 = getStructurePos(type, g->mc, g->seed, regX, regZ, &mpos1);
	if(!monum1)
		return false;
	
	{
		int monum2 = getStructurePos(type, g->mc, g->seed, regX + 1, regZ, &mpos2);
		if(!monum2)
			return false;
		
		if(distpos(mpos1, mpos2) > 2*radius)
			return false;
		
		int viable = isViableStructurePos(type, g, mpos1.x, mpos1.z, 0)
						&& isViableStructurePos(type, g, mpos2.x, mpos2.z, 0);
		if(!viable)
			return false;
	}
	
	{
		int monum2 = getStructurePos(type, g->mc, g->seed, regX, regZ + 1, &mpos2);
		if(!monum2)
			return false;
		
		if(distpos(mpos1, mpos2) > 2*radius)
			return false;
		
		int viable = isViableStructurePos(type, g, mpos1.x, mpos1.z, 0)
						&& isViableStructurePos(type, g, mpos2.x, mpos2.z, 0);
		if(!viable)
			return false;
	}
	
	*oppos = (Pos){ (mpos1.x + mpos2.x)/2, (mpos1.z + mpos2.z)/2 };
	return true;
}

static bool hasDoubleStructure(StructureConfig *config, const int DIST, int radius, Generator *g, Pos pos, Pos *oppos)
{
	int span = 16 * config->regionSize;

	int rdist = DIST / span;
	int rx = pos.x / span;
	int rz = pos.z / span;
	
	int x = 0, z = 0, dirx = 0, dirz = -1;
	for(int c = 0; c < (2 * rdist + 1) * (2 * rdist + 1); ++c)
	{
		if(hasRegDoubleStructure(config->structType, g, rx + x, rz + z, radius, oppos))
			return true;
		
		x += dirx;
		z += dirz;
		if(dirz == -1 ? (x == z + 1) : (abs(x) == abs(z)))
		{
			int ndirx = -dirz;
			dirz = dirx;
			dirx = ndirx;
		}
	}

	return false;
}

static bool hasRegQuadStructure(enum StructureType type, Generator *g, int regX, int regZ, int radius, Pos *oppos)
{
	Pos mpos1, mpos2, mpos3, mpos4;
	
	int monum1 = getStructurePos(type, g->mc, g->seed, regX, regZ, &mpos1);
	if(!monum1)
		return false;
	
	int monum2 = getStructurePos(type, g->mc, g->seed, regX + 1, regZ, &mpos2);
	if(!monum2)
		return false;
	
	int monum3 = getStructurePos(type, g->mc, g->seed, regX, regZ + 1, &mpos3);
	if(!monum3)
		return false;
	
	int monum4 = getStructurePos(type, g->mc, g->seed, regX + 1, regZ + 1, &mpos4);
	if(!monum4)
		return false;
	
	if(!haveMaxEnclosure(mpos1, mpos2, mpos3, mpos4, radius))
		return false;
	
	int viable = isViableStructurePos(type, g, mpos1.x, mpos1.z, 0)
					&& isViableStructurePos(type, g, mpos2.x, mpos2.z, 0)
					&& isViableStructurePos(type, g, mpos3.x, mpos3.z, 0)
					&& isViableStructurePos(type, g, mpos4.x, mpos4.z, 0);
	if(viable)
	{
		*oppos = mpos1;
		return true;
	}
	
	return false;
}

static bool hasQuadStructure(StructureConfig *config, const int DIST, int radius, Generator *g, Pos pos, Pos *oppos)
{
	int span = 16 * config->regionSize;

	int rdist = DIST / span;
	int rx = pos.x / span;
	int rz = pos.z / span;
	
	int x = 0, z = 0, dirx = 0, dirz = -1;
	for(int c = 0; c < (2 * rdist + 1) * (2 * rdist + 1); ++c)
	{
		if(hasRegQuadStructure(config->structType, g, rx + x, rz + z, radius, oppos))
			return true;
		
		x += dirx;
		z += dirz;
		if(dirz == -1 ? (x == z + 1) : (abs(x) == abs(z)))
		{
			int ndirx = -dirz;
			dirz = dirx;
			dirx = ndirx;
		}
	}

	return false;
}

static int dist(int x1, int z1, int x2, int z2)
{
	return (int)hypot(x1 - x2, z1 - z2);
}

static int distpos(Pos pos1, Pos pos2)
{
	return dist(pos1.x, pos1.z, pos2.x, pos2.z);
}

typedef struct
{
	uint64_t seedstart;
	uint64_t seedincr;
	FILE *csv;
} seed_search_data_t;

static void *searchSeeds(void *dataptr)
{
	seed_search_data_t *data = (seed_search_data_t *)dataptr;
	uint64_t seedstart = data->seedstart;
	uint64_t seedincr = data->seedincr;
	FILE *csv = data->csv;
	
	Generator g, ng;
	setupGenerator(&g, MC, 0);
	setupGenerator(&ng, MC, 0);
	
	StrongholdIter sh;
	Pos found_pos, nether_pos, outp_pos;
	Pos biomes_pos[6];
	Pos slimes_pos;
	Pos huts_pos;
	Pos monums_pos;
	int slimecount;
	
	for(uint64_t seed = seedstart;; seed += seedincr)
	{
		if(seed % 10000 == 0)
			printf("Searching %" PRId64 "...\n", seed);
		
		initFirstStronghold(&sh, MC, seed);

		applySeed(&g, 0, seed);
		// applySeed(&ng, -1, seed);
		bool init_nether = false;

		bool found = false;
		for(int i = 1; i <= STRONGHOLDS && !found; ++i)
		{
			if (nextStronghold(&sh, &g) <= 0)
				break;

			// test whether under mushroom island
			if(!isMushroom(&g, sh.pos)) continue;
			if(!init_nether) applySeed(&ng, -1, seed);
			if(!hasFortress(&ng, sh.pos.x / 8, sh.pos.z / 8, &nether_pos)) continue;
			if(!hasOutpost(&g, sh.pos, &outp_pos)) continue;
			if(!hasBiomes(&g, sh.pos, biomes_pos)) continue;
			if((slimecount = findSlimes(seed, sh.pos, &slimes_pos)) == 0) continue;
			if(!(
				count_hut == 2
				? hasDoubleStructure(&conf_hut, HUTS_DIST, 128, &g, sh.pos, &huts_pos)
				: hasQuadStructure(&conf_hut, HUTS_DIST, 128, &g, sh.pos, &huts_pos)
			)) continue;
			if(!(
				count_monum == 2
				? hasDoubleStructure(&conf_monument, MONUMENTS_DIST, 100, &g, sh.pos, &monums_pos)
				: hasQuadStructure(&conf_monument, MONUMENTS_DIST, 100, &g, sh.pos, &monums_pos)
			)) continue;
			
			found = true;
			found_pos = sh.pos;
		}

		if(found)
		{
			printf("Found seed %" PRId64 "\n", seed);
			fflush(stdout);
			// this isn't really thread-safe, but with how rarely seeds are found, that probably doesn't matter
			// (could replace it by multiple sprintf and one fprintf if it _did_ matter)
			fprintf(csv, "%" PRId64, seed);
			fprintf(csv, ";%d;%d;%d", found_pos.x, found_pos.z, countMushs(&g, found_pos));
			fprintf(csv, ";%d;%d;%d", nether_pos.x, nether_pos.z, dist(found_pos.x / 8, found_pos.z / 8, nether_pos.x, nether_pos.z));
			fprintf(csv, ";%d;%d;%d", outp_pos.x, outp_pos.z, distpos(found_pos, outp_pos));
			fprintf(csv, ";%d;%d", biomes_pos[0].x, biomes_pos[0].z);
			fprintf(csv, ";%d;%d", biomes_pos[1].x, biomes_pos[1].z);
			fprintf(csv, ";%d;%d", biomes_pos[2].x, biomes_pos[2].z);
			fprintf(csv, ";%d;%d", biomes_pos[5].x, biomes_pos[5].z);
			fprintf(csv, ";%d;%d;%d", slimes_pos.x, slimes_pos.z, slimecount);
			fprintf(csv, ";%d;%d", huts_pos.x, huts_pos.z);
			fprintf(csv, ";%d;%d", monums_pos.x, monums_pos.z);
			fprintf(csv, "\n");
			fflush(csv);
		}
	}
	
	return NULL;
}

/*
DONE:
- A stronghold on top of a mushroom biome 
- All distances are measured relative to this.
- mushroom size of selected (number of blocks that are mushroom biome within a 512r circle) 
- A nether fortress inside a soul sand valley biome within 1000 blocks (as measured in the nether) of the hub.
- An outpost within 10k of the hub.
- Large enough jungle (within X of hub) (1k blocks by 1k blocks of any jungle type)
- Large enough taiga (within X of hub) (1k by 1k)
- 50+ slime chunks within a 128r circle (despawn sphere) of a point with that point being within 10k of base (nether blocks) - not sure if possible
- Runs continuously finding as many as possible
- x and z coordinates of all things found
- x and z of the center point for the slime chunk as well as the number of slime chunks
- Quad hut (first within X of hub)
- Quad monument (within X of hub)

NOT DONE:
- Mesa region suitable for a quarry (very long, length being min 2000, and at least 5x width) (within X of hub) being at least Y% mesa biome variant
  -> (is 2000x10000 possible???)
  -> currently testing 2000x8
  -> not currently counting specific variant

NOT DOABLE:
- 5-6 crossroad fortress
- First generated end pillar (east) is the widest one - not sure if possible
*/

int main(int argc, char *argv[])
{
	if(argc >= 2 && strcmp(argv[1], "?") == 0)
	{
		printf("Usage: %s [start] [threads] [huts] [monuments]\n", argv[0]);
		printf("       %s ?\n", argv[0]);
		printf("Options:\n");
		printf("  start = seed to start search at, default 0\n");
		printf("  threads = number of threads to use, default %d\n", THREADS);
		printf("  huts = number of witch huts to look for, one of 2,4, default 4\n");
		printf("  monuments = number of monuments to look for, one of 2,4, default 4\n");
		return 0;
	}
	
	uint64_t seedstart = 0;
	if(argc >= 2)
		sscanf(argv[1], "%" PRId64, &seedstart);
	
	int threadcount = THREADS;
	if(argc >= 3)
		sscanf(argv[2], "%d", &threadcount);
	
	if(threadcount <= 0 || threadcount > 64)
	{
		printf("%d is not a good number of threads to use.\n", threadcount);
		return 1;
	}
	
	count_hut = 4;
	if(argc >= 4)
		sscanf(argv[3], "%d", &count_hut);
	
	count_monum = 4;
	if(argc >= 5)
		sscanf(argv[4], "%d", &count_monum);
	
	getStructureConfig(Fortress, MC, &conf_fortress);
	getStructureConfig(Outpost, MC, &conf_outpost);
	getStructureConfig(Swamp_Hut, MC, &conf_hut);
	getStructureConfig(Monument, MC, &conf_monument);

	FILE *csv = fopen("seeds.csv", "a");
	if(!csv)
	{
		printf("Could not open output file.");
		return 1;
	}
	
	char *label[5] = {NULL, NULL, "Double", NULL, "Quad"};
	fprintf(csv, "Seed;Base X;Base Z;Mushroom Count;Fortress X;Fortress Z;Fortress Dist;Outpost X;Outpost Z;Outpost Dist;"
		"Jungle X;Jungle Z;Mega Taiga X;Mega Taiga Z;Mesa Quarry X;Mesa Quarry Z;"
		"Desert X;Desert Z;Slime X;Slime Z;Slime Count;"
		"%s Hut X;%s Hut Z;%s Monument X;%s Monument Z\n",
		label[count_hut], label[count_hut], label[count_monum], label[count_monum]);

	pthread_t threads[64];
	seed_search_data_t datas[64];
	int err = 0;
	for(int i = 0; i < threadcount; ++i)
	{
		datas[i] = (seed_search_data_t) { seedstart + (uint64_t)i, (uint64_t)threadcount, csv };
		err |= pthread_create(&threads[i], NULL, &searchSeeds, &datas[i]);
	}
	if(err)
	{
		printf("Thread creation failed.\n");
		fclose(csv);
		return 1;
	}

	for(int i = 0; i < threadcount; ++i)
	{
		err = pthread_join(threads[i], NULL);
		if(err)
			printf("Thread %d failed.\n", i);
	}
	
	printf("Done!\n");
	fclose(csv);

	return 0;
}
