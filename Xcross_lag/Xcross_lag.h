// Xcross_lag.h : PROJECT_NAME ���ε{�����D�n���Y��
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�� PCH �]�t���ɮ׫e���]�t 'stdafx.h'"
#endif

#include "resource.h"		// �D�n�Ÿ�


// CXcross_lagApp:
// �аѾ\��@�����O�� Xcross_lag.cpp
//

class CXcross_lagApp : public CWinApp
{
public:
	CXcross_lagApp();

// �мg
	public:
	virtual BOOL InitInstance();

// �{���X��@

	DECLARE_MESSAGE_MAP()
};

extern CXcross_lagApp theApp;