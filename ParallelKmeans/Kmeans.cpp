#define _CRT_SECURE_NO_DEPRECATE
#include "PointsCalc.h"
#include "Kmeans.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static void initiateClusters(Point *points, Cluster *clusters, long numPoints, long numClusters);
static void resetClusters(Cluster *clusters, long numClusters);

KmeansAns* runKmeans(Point *points, long numPoints, long numClusters, long maxIter, double step)
{
	KmeansAns *ans = (KmeansAns *)malloc(sizeof(KmeansAns));
	long clusterIndex;
	double newCenterX = 0, newCenterY = 0;
	long i, j;
	bool changed = true;
	double generalMinDist = -1;
	long iter = 0;

	Cluster *clusters = (Cluster *)malloc(numClusters * sizeof(Cluster));

	initiateClusters(points, clusters, numPoints, numClusters);

	while (changed && iter < maxIter) {

		if (iter != 0) {
			resetClusters(clusters, numClusters);
		}

		// calcDistances();
		for (i = 0; i < numPoints; i++) {
			double d = 0, minD = -1;
			for (j = 0; j < numClusters; j++) {
				d = fabs(sqrt(pow(points[i].x - clusters[j].center.x, 2) + pow(points[i].y - clusters[j].center.y, 2)));
				if (minD == -1 || d < minD) {
					minD = d;
					clusterIndex = j;
				}
			}
			//assign point to cluster
			clusters[clusterIndex].clustPoints[clusters[clusterIndex].numClustPoints] = points[i];
			clusters[clusterIndex].numClustPoints++;
		}

		// calcCenters();
		changed = false;
		
		#pragma omp parallel default(none) shared(clusters, changed, numClusters) private(i, j, newCenterX, newCenterY) 
		{
			bool isChanged = false;
			#pragma omp for 
			for (i = 0; i < numClusters; i++) {
				newCenterX = 0;
				newCenterY = 0;
				for (j = 0; j < clusters[i].numClustPoints; j++) {
					newCenterX += clusters[i].clustPoints[j].x;
					newCenterY += clusters[i].clustPoints[j].y;
				}
				newCenterX = newCenterX / clusters[i].numClustPoints;
				newCenterY = newCenterY / clusters[i].numClustPoints;
				if ((clusters[i].center.x != newCenterX) && (clusters[i].center.y != newCenterY)) {
					clusters[i].center.x = newCenterX;
					clusters[i].center.y = newCenterY;
					isChanged = true;
				}
			}
			// no need for critical section since changed is either updated to true or not
			if(isChanged) changed = true;
		}
		iter++;
	}
	// checkClusterMinDist() - cannot parallel with openmp since update of shard var is done in inner loor - a synch will make it sequential and add openmp overhead
	generalMinDist = -1;
	for (i = 0; i < numClusters; i++) {
		double d;
		for (j = i + 1; j < numClusters; j++) {
			d = fabs(sqrt(pow(clusters[i].center.x - clusters[j].center.x, 2) + pow(clusters[i].center.y - clusters[j].center.y, 2)));
			if (generalMinDist == -1 || d < generalMinDist) {
				generalMinDist = d;
			}
		}
	}

	// Kmeans Answer
	(*ans).minDistance = generalMinDist;
	(*ans).timeStep = step;
	(*ans).CentersX = (double*)malloc(numClusters * sizeof(double));
	(*ans).CentersY = (double*)malloc(numClusters * sizeof(double));
	for (i = 0; i < numClusters; i++) {
		(*ans).CentersX[i] = clusters[i].center.x;
		(*ans).CentersY[i] = clusters[i].center.y;
	}

	return ans;
}

static void initiateClusters(Point *points, Cluster *clusters, long numPoints, long numClusters)
{
	long i = 0;

	for (i = 0; i < numClusters; i++) {
		clusters[i].center = points[i];
		clusters[i].index = i;
		clusters[i].numClustPoints = 0;
		clusters[i].clustPoints = (Point *)malloc(numPoints * sizeof(Point));
	}
}

static void resetClusters(Cluster *clusters, long numClusters)
{
	long i = 0;

	for (i = 0; i < numClusters; i++) {
		clusters[i].numClustPoints = 0;
	}
}