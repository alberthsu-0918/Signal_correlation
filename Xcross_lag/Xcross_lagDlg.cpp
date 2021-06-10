// Xcross_lagDlg.cpp : 實作檔
//
#include <windows.h>
#include "stdafx.h"
#include "Xcross_lag.h"
#include "Xcross_lagDlg.h"

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include "fft.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 對 App About 使用 CAboutDlg 對話方塊

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// 對話方塊資料
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支援

// 程式碼實作
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CXcross_lagDlg 對話方塊




CXcross_lagDlg::CXcross_lagDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CXcross_lagDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CXcross_lagDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CXcross_lagDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDOK, &CXcross_lagDlg::OnBnClickedOk)
	ON_BN_CLICKED(btn_fastCCR, &CXcross_lagDlg::OnBnClickedfastccr)
END_MESSAGE_MAP()


// CXcross_lagDlg 訊息處理常式

BOOL CXcross_lagDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// 將 [關於...] 功能表加入系統功能表。

	// IDM_ABOUTBOX 必須在系統命令範圍之中。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 設定此對話方塊的圖示。當應用程式的主視窗不是對話方塊時，
	// 框架會自動從事此作業
	SetIcon(m_hIcon, TRUE);			// 設定大圖示
	SetIcon(m_hIcon, FALSE);		// 設定小圖示

	// TODO: 在此加入額外的初始設定

	return TRUE;  // 傳回 TRUE，除非您對控制項設定焦點
}

void CXcross_lagDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// 如果將最小化按鈕加入您的對話方塊，您需要下列的程式碼，
// 以便繪製圖示。對於使用文件/檢視模式的 MFC 應用程式，
// 框架會自動完成此作業。

void CXcross_lagDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 繪製的裝置內容

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 將圖示置中於用戶端矩形
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 描繪圖示
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// 當使用者拖曳最小化視窗時，
// 系統呼叫這個功能取得游標顯示。
HCURSOR CXcross_lagDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

float XCorr (float *x, float *y, int n, int delay) //n = data number
{
	int i, k;
	float sum = 0;
	float mean1 = 0, mean2 = 0;

	/*
	for (i = 0; i < n; i++)
	{
		mean1 += x[i];
		mean2 += y[i];
	}
	mean1 /= n;
	mean2 /= n;
	*/

	for (i = 0; i < n; i++)
	{
		k = i + delay;
		if(k >= n)
			continue;
		else
			sum += (x[i]) * (y[k]);
		//k = (i + delay) % n;
		//sum += (x[i] - mean1) * (y[k] - mean2);
		//sum += (x[i]) * (y[k]);
	}

	return (sum);
}

float* Fast_XCorr (float *x, float *y, unsigned long n, float *ZR, float *ZI) //n = data number
{
	int k;
	float sum = 0;
	float mean1 = 0, mean2 = 0;
	unsigned int err = 0;
	float *X_real, *X_img, *Y_real, *Y_img;
	float *z_R, *z_I;

	z_R = new float[n];
	z_I = new float[n];

	X_real = new float[n];
	X_img = new float[n];

	err = fft_float(
      0   // [IN ] 0) forward, x) inverse transform
    , (unsigned int)n   // [IN ] number of samples, must be power of 2
    , x   // [IN ] source samples, real part
    , 0   // [IN ] source samples, image part, could be 0
    , X_real   // [OUT] target samples, real part
    , X_img   // [OUT] target samples, image part
    );

	Y_real = new float[n];
	Y_img = new float[n];

	err = fft_float(
      0   // [IN ] 0) forward, x) inverse transform
    , (unsigned int)n   // [IN ] number of samples, must be power of 2
    , y   // [IN ] source samples, real part
    , 0   // [IN ] source samples, image part, could be 0
    , Y_real   // [OUT] target samples, real part
    , Y_img   // [OUT] target samples, image part
    );

	for(int i=0;i<n;i++)
	{
		X_img[i] = -1 * X_img[i]; //take conj
		z_R[i] = X_real[i]* Y_real[i] - X_img[i]* Y_img[i];
		z_I[i] = X_img[i]* Y_real[i] +  X_real[i]* Y_img[i];
	}

	delete [] X_real;
	delete [] X_img;
	delete [] Y_real;
	delete [] Y_img;

	//take inverse fft
	err = fft_float( 
      1   // [IN ] 0) forward, x) inverse transform
    , (unsigned int)n   // [IN ] number of samples, must be power of 2
    , z_R   // [IN ] source samples, real part
    , z_I   // [IN ] source samples, image part, could be 0
    , ZR   // [OUT] target samples, real part
    , ZI   // [OUT] target samples, image part
    );

	return (ZR);
}

unsigned long FindMinPowerof2(unsigned long x)
{
	double temp = (double)(unsigned long)(log((F64)x)/log(2.0)); 
	return (unsigned long)(pow(2.0,temp)); 
}

int FindMaxIndex(float *correlation, int length)
{
	int index = 0;
	float ref = 0.0;
	for(int i=0; i<length; i++)
	{
		if(abs(correlation[i]) > ref){
			ref = correlation[i];
			index = i;
		}
	}

	return index;
}

void CXcross_lagDlg::OnBnClickedOk()
{
	// TODO: 在此加入控制項告知處理常式程式碼
	FILE *fp_file1, *fp_file2, *fp_xcfile;
	fp_file1 = fopen("x.txt", "r");
	fp_file2 = fopen("y.txt", "r");
	fp_xcfile = fopen("xc.txt","w");
	int length1 = 0, length2 = 0;
	int j, delayRange;
	float value_ref;
	float *x, *y, *corr;
	int Actual_delay = 0;

	//check file's size
	while((fscanf(fp_file1,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		length1++;
	}
	rewind(fp_file1);

	while((fscanf(fp_file2,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		length2++;
	}
	rewind(fp_file2);

	if(length1!=length2){
		//fuck you
	}

	x = new float[length1];
	y = new float[length2];

	j = 0;
	while((fscanf(fp_file1,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		x[j] = value_ref;
		j++;
	}

	j = 0;
	while((fscanf(fp_file2,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		y[j] = value_ref;
		j++;
	}
/*
	for(int i =0;i<length1;i++){
		x[i+length1] = 0;
		y[i+length1] = 0;
	}
*/
	corr = new float[2*length1];

	delayRange = 65536;

	for(int delay =0; delay< delayRange; delay++)
	{
		corr[delay] = XCorr(x, y, length1, delay); //corr[]: correlation coefficients
		fprintf(fp_xcfile, "%f\n", corr[delay]);
	}

	Actual_delay = FindMaxIndex(corr, 2*length1);

	//OnOK();

	delete [] x;
	delete [] y;
	delete [] corr;
	fclose(fp_file1);
	fclose(fp_file2);
	fclose(fp_xcfile);
}



void CXcross_lagDlg::OnBnClickedfastccr()
{
	// TODO: 在此加入控制項告知處理常式程式碼
	FILE *fp_file1, *fp_file2, *fp_fxcfile;
	//fp_file1 = fopen("x.txt", "r"); //original signal, so called your output or golden pattern
	//fp_file2 = fopen("y.txt", "r"); //acquried sginal, the lag signal
	fp_file1 = fopen("origin.csv", "r"); //original signal, so called your output or golden pattern
	fp_file2 = fopen("rec.csv", "r"); //acquried sginal, the lag signal
	fp_fxcfile = fopen("fast_xc.txt","w");
	unsigned long length1 = 0, length2 = 0, compensationLength = 0;
	int j;
	float value_ref;
	float *x, *y, *corr, *xy_real, *xy_img;
	int Actual_delay = 0;
	CString delayMesg;

	//check file's size
	while((fscanf(fp_file1,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		length1++;
	}
	rewind(fp_file1);

	while((fscanf(fp_file2,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		length2++;
	}
	rewind(fp_file2);

	if(length1<=length2){
		//fuck you
		compensationLength = 2*FindMinPowerof2(length2);
	}
	else
		compensationLength = 2*FindMinPowerof2(length1);


	x = new float[compensationLength];
	memset(x, (int) 0, compensationLength*sizeof(float));
	y = new float[compensationLength];
	memset(y, (int) 0, compensationLength*sizeof(float));

	j = 0;
	while((fscanf(fp_file1,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		x[j] = value_ref;
		j++;
	}

	j = 0;
	while((fscanf(fp_file2,"%f", &value_ref)) !=EOF){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		y[j] = value_ref;
		j++;
	}

	xy_real = new float[compensationLength];
	xy_img = new float[compensationLength];

	corr =  Fast_XCorr (x, y, compensationLength, xy_real, xy_img ); //n = data number, x: original signal, y: lag signal 

	//for debug, save all correlation coefficent
	
	for(int plotNum =0; plotNum < compensationLength; plotNum++)
	{
		fprintf(fp_fxcfile, "%f\n", corr[plotNum]);
	}
	

	Actual_delay = FindMaxIndex(corr, compensationLength);
	delayMesg.Format(_T("Delay betwwen x & y: %d Clocks"), Actual_delay);
	AfxMessageBox(delayMesg);


	delete [] xy_real;
	delete [] xy_img;
	delete [] x;
	delete [] y;
	
	fclose(fp_file1);
	fclose(fp_file2);
	fclose(fp_fxcfile);
}
