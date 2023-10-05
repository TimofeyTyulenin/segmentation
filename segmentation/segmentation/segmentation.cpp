// clustering.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;
int dp = 0, fp = 0;
struct Point {
	double x, y, z;
	double centrx, centry, centrz;
	bool entered = false;
};

//высчитывает расстояние между двумя точками
double distance(const Point& p1, const Point& p2) {
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}

//принимает на вход облако точек, количество итераций и пороговое значение для межквартильного расстояния
void filter(vector<Point>& points, int iterations, double threshold) {
	for (int iter = 0; iter < iterations; ++iter) {
		vector<double> distances;
		//далее находим минимальное расстояние между точками
		for (int i = 0; i < points.size(); ++i) {
			double min_dist = numeric_limits<double>::max(); //присваивается наибольшее возможное	значение для double

			for (int j = 0; j < points.size(); ++j) {
				if (i != j) {
					double dist = distance(points[i], points[j]); //находим расстояние между двумя точками
					//если расстояние между двумя точками меньше выбранного минимального расстояния, то 
					//перезаписываем min_dist
					if (dist < min_dist) {
						min_dist = dist;
					}
				}
			}
			distances.push_back(min_dist); //добавляем минимальную дистанцию от данной точки до другой в вектор значений
		}
		//в итоге получаем вектор из минимальных расстояний до других точек для каждой точки

		sort(distances.begin(), distances.end()); //сортировка значений вектора
		double median = distances[distances.size() / 2]; //находим середину
		double q1 = distances[distances.size() / 4];
		double q3 = distances[3 * distances.size() / 4];
		double iqr = q3 - q1;
		double threshold_distance = median + threshold * iqr;
		vector<Point> filtered_points;

		//фильтруем набор значений
		for (int i = 0; i < points.size(); ++i) {
			double min_dist = numeric_limits<double>::max();
			for (int j = 0; j < points.size(); ++j) {
				if (i != j) {
					double dist = distance(points[i], points[j]);
					if (dist < min_dist) {
						min_dist = dist;
					}
				}
			}
			if (min_dist < threshold_distance) //проверка, что min_dist меньше порогового расстояния
			{
				filtered_points.push_back(points[i]);
				fp++;
			}
			else {
				dp++;
			}
		}
		points = filtered_points;
	}
}

double euclidean_distance(Point p1, Point p2) {
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

vector<Point> mean_shift(vector<Point>& points, double radius) {

	vector<Point> centroids;
	vector<bool> visited(points.size(), false);
	vector<Point> clustered;
	Point new_centroid;
	double e = radius; //точность
	const int ps = points.size();
	int* mas = new int[ps];
	int l = 0;


	for (int y = 0; y < points.size(); y++)
	{
		for (int g = 0; g < points.size() - y - 1; g++)
		{
			if (points[g].z > points[g + 1].z || points[g].x > points[g + 1].x); 
			{
				swap(points[g], points[g + 1]);
			}
		}

	}

	


	for (int i = 0; i < points.size(); i++) {
		if (visited[i]) continue; //если посещали эту точку, то пропускаем итерацию


		Point centroid = points[i];
		int count = 0;

		while (true) {
			Point sum = { 0, 0, 0 };

			//поиск среднего в области "подъем на гору"
			for (int j = 0; j < points.size(); j++) {
				if (euclidean_distance(points[j], centroid) <= radius) {
					sum.x += points[j].x;
					sum.y += points[j].y;
					sum.z += points[j].z;


					if (!visited[j]) {
						mas[l] = j;
						l++;
						clustered.push_back(points[j]);
					}
					visited[j] = true;



					count++;
				}
			}

			if (sum.x != 0 && sum.y != 0 && sum.z != 0)
				new_centroid = { sum.x / count, sum.y / count, sum.z / count }; //расчет новой центроиды

			//момент, когда пора остановиться, выдохнуть и попить пивка
			if (euclidean_distance(new_centroid, centroid) < e) {
				centroids.push_back(new_centroid);
				for (int a = 0; a < clustered.size(); a++)
				{
					int n = mas[a];
					points[n].centrx = new_centroid.x;
					points[n].centry = new_centroid.y;
					points[n].centrz = new_centroid.z;

				}
				for (int a = 0; a < clustered.size(); a++)
					mas[a] = 0;
				l = 0;
				clustered.clear();


				break;
			}


			centroid = new_centroid;

		}
	}

	return centroids;
}









int main() {
	string file = "lep3.txt";
	ifstream fin;
	Point** pfinally;
	vector<Point> points = { };
	int count = 0, ncl = 0, centrsize = 0;
	fin.open(file);
	if (!fin.is_open())
	{
		cout << "Error" << endl;
		return 9;
	}
	else {
		cout << "File opened" << endl;

		while (!fin.eof()) {

			double x, y, z;
			fin >> x >> y >> z;
			points.push_back({ x,y,z });
			count++;
		}
		fin.close();

	}
	filter(points, 1, 33);
	cout << "dp = " << dp << endl << "count = " << count << endl;

	double radius = 10;

	vector<Point> centroids = mean_shift(points, radius);

	double epsx = 32;
	double epsy = epsx*4;// *7;
	centrsize = centroids.size();

	int* clusters = new int[centroids.size()];
	int cluster = 0;
	int k = 0, t = 0;



	



	int cns = centrsize;


	

	for (int i = 0; i < centrsize; i++)
	{
		for (int j = i + 1; j < centrsize; j++)
		{
			if ((abs(centroids[i].x - centroids[j].x) < epsx) && (abs(centroids[i].y - centroids[j].y) < epsy))
			{
				centroids[j].x = centroids[i].x;
				centroids[j].y = centroids[i].y;
				centroids[j].z = centroids[i].z;

			}
		}
	}
	for (int t = 0; t < 3; t++) {
	for (int j = 0; j < centroids.size(); j++) {
		for (int i = j + 1; i < centroids.size(); i++)
		{
			if (centroids[j].x == centroids[i].x) {
				swap(centroids[i], centroids[centroids.size() - 1]);
				centroids.pop_back();
			}
		}
	}


}
	cns = centroids.size();

	pfinally = new Point * [cns];
	for (int count = 0; count < cns; count++)
		pfinally[count] = new Point[points.size()];







	ncl = 0;
	for (int i = 0; i < cns;i++) {
		int c = 0;

		for (int j = 0; j < points.size(); j++)
		{
			if (((abs(points[j].centrx - centroids[i].x) < epsx) && ((abs(points[j].centry - centroids[i].y) < epsy))) || (points[j].centrz == centroids[i].z))
			{
				if (!points[j].entered)
				{
					pfinally[ncl][c] = points[j];
					points[j].entered = true;
					c++;
				}
				
			}

		}
		ncl++;





		cout << "Centroid: (" << centroids[i].x << ", " << centroids[i].y << ", " << centroids[i].z << ")" << endl;

	}

	

	


	ofstream out;          // поток для записи
	out.open("outlep3.txt");      // открываем файл для записи
	if (out.is_open())
	{
		for (auto centroid : centroids) {
			out << centroid.x << '\t' << centroid.y << '\t' << centroid.z << endl;
		}
		out << "s" << endl;

		for (int v = 0; v < ncl; v++)
		{
			for (int w = 0; w < points.size(); w++)
			{
				if (pfinally[v][w].x > 0) {
					

					out << pfinally[v][w].x << '\t' << pfinally[v][w].y << '\t' << pfinally[v][w].z;
					out << endl;
				}
			}
			out << "s" << endl;
		}

	}
	else {
		cout << "File wasnt written" << endl;
		return 8;
	}
	out.close();
	cout << "File has been written" << endl;

	return 0;
}
