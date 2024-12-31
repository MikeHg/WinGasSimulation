// WinGasSimulation.cpp : Definiert den Einstiegspunkt für die Anwendung.
//

#include "pch.h"
#include "framework.h"
#include "WinGasSimulation.h"
#include "Direct2DApp.h"

#define MAX_LOADSTRING 100

// Globale Variablen:
HINSTANCE hInst;                                // Aktuelle Instanz
WCHAR szTitle[MAX_LOADSTRING];                  // Titelleistentext
WCHAR szWindowClass[MAX_LOADSTRING];            // Der Klassenname des Hauptfensters.

// Vorwärtsdeklarationen der in diesem Codemodul enthaltenen Funktionen:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int, HWND*);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	OptionsH(HWND, UINT, WPARAM, LPARAM);
DWORD WINAPI		ThreadProc(LPVOID);
HWND hwndOptions = NULL;


// Global instance of the Direct2DApp class
Direct2DApp g_Direct2DApp;
HANDLE thrdHndl = NULL;
DWORD threadId = 0;
int volatile stop_thrd = 1;

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Hier Code einfügen.
    box = new GasBox();

    // Globale Zeichenfolgen initialisieren
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_WINGASSIMULATION, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    HWND myhWnd;
    // Anwendungsinitialisierung ausführen:
    if (!InitInstance (hInstance, nCmdShow, &myhWnd))
    {
        return FALSE;
    }

    // Send WM_TIMER event every 40ms so that OnPaint is regularly called
    SetTimer(myhWnd, 4711, 40, NULL);


    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_WINGASSIMULATION));

    //DialogBox(hInstance, MAKEINTRESOURCE(IDD_DIALOG1), myhWnd, OptionsH);
    hwndOptions = CreateDialog(hInstance, MAKEINTRESOURCE(IDD_DIALOG1), myhWnd, OptionsH);
    ShowWindow(hwndOptions, SW_SHOW);

    MSG msg;

    // Hauptnachrichtenschleife:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    KillTimer(myhWnd, 4711);

    return (int) msg.wParam;
}



//
//  FUNKTION: MyRegisterClass()
//
//  ZWECK: Registriert die Fensterklasse.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_WINGASSIMULATION));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_WINGASSIMULATION);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNKTION: InitInstance(HINSTANCE, int)
//
//   ZWECK: Speichert das Instanzenhandle und erstellt das Hauptfenster.
//
//   KOMMENTARE:
//
//        In dieser Funktion wird das Instanzenhandle in einer globalen Variablen gespeichert, und das
//        Hauptprogrammfenster wird erstellt und angezeigt.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow, HWND* retWndH)
{
   hInst = hInstance; // Instanzenhandle in der globalen Variablen speichern

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
       100, 100, (int)(2.0 * BOX_X) + 700 + 1000, (int)(2.0 * BOX_Y) + 610, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   *retWndH = hWnd;

   return TRUE;
}

//
//  FUNKTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  ZWECK: Verarbeitet Meldungen für das Hauptfenster.
//
//  WM_COMMAND  - Verarbeiten des Anwendungsmenüs
//  WM_PAINT    - Darstellen des Hauptfensters
//  WM_DESTROY  - Ausgeben einer Beendenmeldung und zurückkehren
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
    case WM_ERASEBKGND:
        return 1;

    case WM_CREATE:
        g_Direct2DApp.Initialize(hWnd);
        return 0;

    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Menüauswahl analysieren:
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case ID_OPTIONS_MDCONTROL:
                if (!IsWindow(hwndOptions)) {
                    hwndOptions = CreateDialog(hInst, MAKEINTRESOURCE(IDD_DIALOG1), hWnd, OptionsH);
                    ShowWindow(hwndOptions, SW_SHOW);
                }
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        g_Direct2DApp.OnPaint(hWnd);
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_TIMER:
        //if (stop_thrd == 0)
        InvalidateRgn(hWnd, NULL, TRUE);
        break;
    case WM_SIZE:
        g_Direct2DApp.OnResize(hWnd);
        //Resize();
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Meldungshandler für Infofeld.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}

// Message handler for options box.
INT_PTR CALLBACK OptionsH(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	DWORD dwWaitResult;

	switch (message)
	{
	case WM_INITDIALOG:
        if (stop_thrd == 0) {
            EnableWindow(GetDlgItem(hDlg, IDC_BUTTON1), FALSE);
            EnableWindow(GetDlgItem(hDlg, IDC_BUTTON3), FALSE);
        }
		return (INT_PTR)TRUE;

    case WM_CLOSE:
	case WM_DESTROY:
        EndDialog(hDlg, LOWORD(wParam));
        hwndOptions = NULL;
		return (INT_PTR)TRUE;

	case WM_COMMAND:

        // Close window button
		if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
            hwndOptions = NULL;
			return (INT_PTR)TRUE;
		}

        // Start Simulation Button
        if (LOWORD(wParam) == IDC_BUTTON1)
        {
            stop_thrd = 0;
            thrdHndl = CreateThread(NULL, 0, ThreadProc, NULL, 0, &threadId);
            if (threadId == NULL)
                ExitProcess(1);

            EnableWindow(GetDlgItem(hDlg, IDC_BUTTON1), FALSE);
            EnableWindow(GetDlgItem(hDlg, IDC_BUTTON3), FALSE);
        }

        //Stop Simulation Button
        if (LOWORD(wParam) == IDC_BUTTON2)
        {
            if (thrdHndl != NULL)
            {
                // Tell thread to stop and wait for stopping
                stop_thrd = 1;

                // Wait for thread to terminate
                dwWaitResult = WaitForSingleObject(thrdHndl, 5000);
                switch (dwWaitResult)
                {
                case WAIT_OBJECT_0:
                    CloseHandle(thrdHndl);
                    thrdHndl = NULL;
                    EnableWindow(GetDlgItem(hwndOptions, IDC_BUTTON1), TRUE);
                    EnableWindow(GetDlgItem(hwndOptions, IDC_BUTTON3), TRUE);
                    break;
                case WAIT_TIMEOUT:
                case WAIT_ABANDONED:
                    break;
                }
            }
        }

        //Switch velocities
        if (LOWORD(wParam) == IDC_BUTTON3) {
            if(stop_thrd == 1)
                box->switchVelocities();
        }

    default:
        return (INT_PTR)FALSE;
		
	}

	return (INT_PTR)FALSE;
}

DWORD WINAPI ThreadProc(LPVOID lpParam)
{

    for (int i = 0; i < 3500*3; i++) {
        box->verlet();
        if (stop_thrd != 0) return 0;
    }

    if (IsWindow(hwndOptions)) {
        EnableWindow(GetDlgItem(hwndOptions, IDC_BUTTON1), TRUE);
        EnableWindow(GetDlgItem(hwndOptions, IDC_BUTTON3), TRUE);
    }

    stop_thrd = 1;
    return 0;
}
