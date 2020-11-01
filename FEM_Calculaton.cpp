#include<iostream>
#include<vector>
#include<fstream>

#define xNodeNumber 41
#define yNodeNumber 41
#define dx 0.01
#define dy 0.01
#define eps0 8.85e-12
#define rho 0
#define potential 100

using namespace std;

int main()
{
	int allNodeNumber = xNodeNumber * yNodeNumber;
	int allElementNumber = (xNodeNumber - 1) * (yNodeNumber - 1) * 2;
	int i, j;
	int nNode = 0;
	int nElement = 0;
	double Area = 0.5 * dx * dy;

	//定义数组
	double Ke[3][3], be[3];
	//double xElementNode[3], yElementNode[3];
	double bn[3], cn[3];
	double* xAllNode = new double[allNodeNumber];
	double* yAllNode = new double[allNodeNumber];
	int** allNodeIdentifier = new int*[allNodeNumber];
	for (i = 0; i < xNodeNumber; i++)
	{
		allNodeIdentifier[i] = new int[yNodeNumber];
	}
	int* allElementIdentifier = new int[allElementNumber];
	int** elementNodeIdentifier = new int* [3];
	for (i = 0; i < 3; i++)
	{
		elementNodeIdentifier[i] = new int[allElementNumber];
	}
	double** K = new double* [allNodeNumber];
	for (i = 0; i < allNodeNumber; i++)
	{
		K[i] = new double[allNodeNumber];
	}
	double* b = new double[allNodeNumber];
	double* U = new double[allNodeNumber];
	double** Aug = new double* [allNodeNumber];
	for (i = 0; i < allNodeNumber; i++)
	{
		Aug[i] = new double[allNodeNumber + 1];
	}

	//数组初始化
	for (i = 0; i < allNodeNumber; i++)
	{
		for (j = 0; j < allNodeNumber; j++)
		{
			K[i][j] = 0;
		}
		b[i] = 0;
	}

	//所有节点编号
	for (j = 0; j < yNodeNumber; j++)
	{
		for (i = 0; i < xNodeNumber; i++)
		{
			allNodeIdentifier[i][j] = nNode;
			xAllNode[nNode] = i * dx;
			yAllNode[nNode] = j * dy;
			nNode++;
		}
	}

	//单元节点编号对应所有节点编号
	for (j = 0; j < yNodeNumber - 1; j++)
	{
		for (i = 0; i < xNodeNumber - 1; i++)
		{
			elementNodeIdentifier[0][nElement] = allNodeIdentifier[i][j + 1];
			elementNodeIdentifier[1][nElement] = allNodeIdentifier[i][j];
			elementNodeIdentifier[2][nElement] = allNodeIdentifier[i + 1][j];
			nElement++;
			elementNodeIdentifier[0][nElement] = allNodeIdentifier[i][j + 1];
			elementNodeIdentifier[1][nElement] = allNodeIdentifier[i + 1][j];
			elementNodeIdentifier[2][nElement] = allNodeIdentifier[i + 1][j + 1];
			nElement++;
		}
	}

	//计算Ke,be,并写入K,b
	for (nElement = 0; nElement < allElementNumber; nElement++)
	{
		int n0, n1, n2;
		//double a, b, c;
		n0 = elementNodeIdentifier[0][nElement];
		n1 = elementNodeIdentifier[1][nElement];
		n2 = elementNodeIdentifier[2][nElement];
		bn[0] = yAllNode[n1] - yAllNode[n2];
		bn[1] = yAllNode[n2] - yAllNode[n0];
		bn[2] = yAllNode[n0] - yAllNode[n1];
		cn[0] = xAllNode[n2] - xAllNode[n1];
		cn[1] = xAllNode[n0] - xAllNode[n2];
		cn[2] = xAllNode[n1] - xAllNode[n0];
		for (i = 0; i < 3; i++)
		{
			int iTemp = elementNodeIdentifier[i][nElement];
			for (j = 0; j < 3; j++)
			{
				Ke[i][j] = (bn[i] * bn[j] + cn[i] * cn[j]) / (4 * Area);
				//cout << Ke[i][j] << endl;
				int jTemp = elementNodeIdentifier[j][nElement];
				K[iTemp][jTemp] += Ke[i][j];
			}
			be[i] = Area / 3 * rho / eps0;
			b[iTemp] += be[i];
		}
	}

	//设置边界条件
	vector<int> p;
	int temp;
	for (i = 0; i < xNodeNumber; i++)
	{
		temp = allNodeIdentifier[i][0];
		p.push_back(temp);
		b[temp] = 0;
		//temp = allNodeIdentifier[i][yNodeNumber - 1];
		//p.push_back(temp);
		//b[temp] = 0;
	}
	for (j = 1; j < xNodeNumber - 1; j++)
	{
		temp = allNodeIdentifier[xNodeNumber - 1][j];
		p.push_back(temp);
		b[temp] = 0;
	}
	for (i = 0; i < int(xNodeNumber / 2); i++)
	{
		temp = allNodeIdentifier[i][int(yNodeNumber / 2)];
		p.push_back(temp);
		b[temp] = potential;
		//temp = allNodeIdentifier[i][int(yNodeNumber / 2) - 1];
		//p.push_back(temp);
		//b[temp] = potential;
	}
	int z = p.size();
	for (i = 0; i < z; i++)
	{
		temp = p[i];
		K[temp][temp] = 1;
		for (j = 0; j < allNodeNumber; j++)
		{
			if (j != temp)
			{
				//b[j] -= K[j][temp] * b[i];
				K[temp][j] = 0;
				//K[j][temp] = 0;
			}
		}
	}
	//for (i = 0; i < allNodeNumber; i++)
	//{
	//	//for (j = 0; j < allNodeNumber; j++)
	//	//{
	//	//	cout << K[i][j] << "   ";
	//	//}
	//	cout << b[i] << endl;
	//}

	//求解KU=b
	int i1, k;
	double max, l, s, temp1;
	for (i = 0; i < allNodeNumber; i++)
	{
		for (j = 0; j < allNodeNumber; j++)          //构建增广矩阵Aug
		{
			Aug[i][j] = K[i][j];
			Aug[i][allNodeNumber] = b[i];
		}
	}
	for (i1 = 0; i1 < allNodeNumber - 1; i1++)                 //求主元
	{
		max = fabs(Aug[i1][i1]);
		k = i1;
		for (i = i1; i < allNodeNumber; i++)
		{
			if (max < fabs(Aug[i][i1]))
			{
				max = fabs(Aug[i][i1]);
				k = i;
			}
		}
		for (j = i1; j < allNodeNumber + 1; j++)       //换行（将最大元素行换到目标行）
		{
			temp1 = Aug[i1][j];
			Aug[i1][j] = Aug[k][j];
			Aug[k][j] = temp1;
		}
		for (k = i1 + 1; k < allNodeNumber; k++)              //消元
		{
			l = -Aug[k][i1] / Aug[i1][i1];
			for (j = i1; j < allNodeNumber + 1; j++)
				Aug[k][j] = Aug[k][j] + l * Aug[i1][j];
		}
	}
	U[allNodeNumber - 1] = Aug[allNodeNumber - 1][allNodeNumber] / Aug[allNodeNumber - 1][allNodeNumber - 1];    //回代求解
	for (i = allNodeNumber - 2; i >= 0; i = i - 1)
	{
		s = 0;
		for (j = i + 1; j < allNodeNumber; j++)
			s = s + Aug[i][j] * U[j];
		U[i] = (Aug[i][allNodeNumber] - s) / Aug[i][i];
		//cout << U[i] << endl;
	}

	//输出结果
	ofstream fp("Potential.dat", ofstream::out);
	fp << "title=\"Potential\"" << endl;
	fp << "variables = \"x\",\"y\",\"U(V)\"" << endl;
	fp << "zone i=" << xNodeNumber << ",j=" << yNodeNumber << endl;
	for (j = 0; j < yNodeNumber; j++)
	{
		for (i = 0; i < xNodeNumber; i++)
		{
			fp << i * dx << " " << j * dx << " " << U[i + j * xNodeNumber] << endl;
		}
	}
	fp.close();
}