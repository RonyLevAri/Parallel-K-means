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

		printf("Iter #%ld\n", iter);

		if (iter != 0) {
			resetClusters(clusters, numClusters);
		}


		generalMinDist = -1;

		// calcDistances();
		for (i = 0; i < numPoints; i++) {
			printf("point #%ld: x = %lf y = %lf\n", i, points[i].x, points[i].y);
			double d = 0, minD = -1;
			for (j = 0; j < numClusters; j++) {
				printf("cluster #%ld with center: x = %lf y = %lf\n", j, clusters[j].center.x, clusters[j].center.y);
				d = fabs(sqrt(pow(points[i].x - clusters[j].center.x, 2) + pow(points[i].y - clusters[j].center.y, 2)));
				printf("calculated distance from center %ld (%lf, %lf) = %lf calculated numD = %lf\n", j, clusters[j].center.x, clusters[j].center.y, d, minD);
				if (minD == -1 || d < minD) {
					minD = d;
					clusterIndex = j;
				}
			}
			printf("Point #%ld is in center #%ld\n", i, clusterIndex);
			//assign point to cluster
			clusters[clusterIndex].clustPoints[clusters[clusterIndex].numClustPoints] = points[i];
			clusters[clusterIndex].numClustPoints++;
			printf("Cluster #%ld has %ld points\n", clusterIndex, clusters[clusterIndex].numClustPoints);
			if ((generalMinDist == -1) || (minD < generalMinDist)) {
				generalMinDist = minD;
			}
			printf("General min distance is %lf\n", generalMinDist);
		}

		// calcCenters();
		changed = false;
		for (i = 0; i < numClusters; i++) {
			newCenterX = 0;
			newCenterY = 0;
			printf("cluster #%ld with center: x = %lf y = %lf\n", i, clusters[i].center.x, clusters[i].center.y);
			for (j = 0; j < clusters[i].numClustPoints; j++) {
				newCenterX += clusters[i].clustPoints[j].x;
				newCenterY += clusters[i].clustPoints[j].y;
			}
			newCenterX = newCenterX / clusters[i].numClustPoints;
			newCenterY = newCenterY / clusters[i].numClustPoints;
			printf("cluster #%ld new center: x = %lf y = %lf\n", i, newCenterX, newCenterY);
			if ((clusters[i].center.x != newCenterX) && (clusters[i].center.y != newCenterY)) {
				changed = true;
				clusters[i].center.x = newCenterX;
				clusters[i].center.y = newCenterY;
			}
		}
		iter++;
	}
	// Kmeans Answer
	(*ans).minDistance = generalMinDist;
	(*ans).timeStep = step;
	(*ans).centers = (Point*)malloc(numClusters * sizeof(Point));
	for (i = 0; i < numClusters; i++) {
		(*ans).centers[i] = clusters[i].center;
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