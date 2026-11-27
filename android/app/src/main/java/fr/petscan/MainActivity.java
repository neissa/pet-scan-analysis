package fr.petscan;

import android.app.Activity;
import android.app.AlertDialog;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.provider.DocumentsContract;
import android.provider.Settings;
import android.view.Gravity;
import android.view.View;
import android.webkit.WebResourceRequest;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.Button;
import android.widget.FrameLayout;

import com.chaquo.python.PyObject;
import com.chaquo.python.Python;
import com.chaquo.python.android.AndroidPlatform;

import org.json.JSONObject;

import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.Socket;

/** Lance le serveur Python (viewer.py) dans l'appli et l'affiche dans une WebView plein écran. */
public class MainActivity extends Activity {
    static final int PORT = 8000;
    static final String URL = "http://127.0.0.1:" + PORT + "/";
    static final int PICK_FOLDER = 1;
    static boolean serverStarted;  // un seul serveur par processus, même si l'activité est recréée

    WebView web;

    @Override
    protected void onCreate(Bundle saved) {
        super.onCreate(saved);
        FrameLayout root = new FrameLayout(this);
        web = new WebView(this);
        WebSettings s = web.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setAllowFileAccess(false);
        s.setAllowContentAccess(false);
        s.setMixedContentMode(WebSettings.MIXED_CONTENT_NEVER_ALLOW);
        web.setWebViewClient(new WebViewClient() {
            @Override
            public boolean shouldOverrideUrlLoading(WebView v, WebResourceRequest r) {
                return !"127.0.0.1".equals(r.getUrl().getHost());  // aucune navigation hors du serveur local
            }
        });
        root.addView(web);
        Button folder = new Button(this);
        folder.setText("Dossier");
        folder.setAlpha(0.75f);
        folder.setOnClickListener(v -> pickFolder());
        root.addView(folder, new FrameLayout.LayoutParams(FrameLayout.LayoutParams.WRAP_CONTENT,
                FrameLayout.LayoutParams.WRAP_CONTENT, Gravity.BOTTOM | Gravity.END));
        setContentView(root);
        web.loadData("<p style='font:20px sans-serif'>Démarrage du serveur…</p>", "text/html; charset=utf-8", null);
        askStorage();
        new Thread(this::startServer).start();
    }

    void startServer() {
        try {
            synchronized (MainActivity.class) {
                if (!serverStarted) {
                    File page = new File(getFilesDir(), "index.html");
                    try (InputStream in = getAssets().open("index.html"); OutputStream out = new FileOutputStream(page)) {
                        byte[] b = new byte[65536];
                        for (int n; (n = in.read(b)) > 0; ) out.write(b, 0, n);
                    }
                    if (!Python.isStarted()) Python.start(new AndroidPlatform(this));
                    Python py = Python.getInstance();
                    PyObject env = py.getModule("os").get("environ");
                    env.callAttr("__setitem__", "PETSCAN_INDEX", page.getAbsolutePath());
                    env.callAttr("__setitem__", "PETSCAN_CONFIG", new File(getFilesDir(), "config.json").getAbsolutePath());
                    py.getModule("sys").put("argv", py.getBuiltins().callAttr("list",
                            (Object) new String[]{"viewer", "--port", String.valueOf(PORT)}));
                    Thread t = new Thread(() -> py.getModule("viewer").callAttr("main"), "serveur-python");
                    t.setDaemon(true);
                    t.start();
                    serverStarted = true;
                }
            }
            for (int i = 0; i < 600; i++) {  // 60 s au plus (premier lancement : extraction des modules)
                try (Socket so = new Socket()) {
                    so.connect(new InetSocketAddress("127.0.0.1", PORT), 200);
                    runOnUiThread(() -> web.loadUrl(URL));
                    return;
                } catch (Exception e) {
                    Thread.sleep(100);
                }
            }
            throw new RuntimeException("le serveur ne répond pas sur le port " + PORT);
        } catch (Exception e) {
            String msg = "<p style='font:18px sans-serif'>Erreur au démarrage : " + android.text.Html.escapeHtml(String.valueOf(e)) + "</p>";
            runOnUiThread(() -> web.loadData(msg, "text/html; charset=utf-8", null));
        }
    }

    /** Le serveur lit les DICOM par leur chemin : il faut l'accès au stockage partagé. */
    void askStorage() {
        if (Build.VERSION.SDK_INT >= 30) {
            if (Environment.isExternalStorageManager()) return;
            new AlertDialog.Builder(this)
                    .setTitle("Accès aux fichiers")
                    .setMessage("Pour lire les dossiers DICOM (Pictures, DCIM, Download…), autorisez PetScan à accéder à tous les fichiers. Les données restent sur la tablette.")
                    .setPositiveButton("Autoriser", (d, w) -> startActivity(new Intent(Settings.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION,
                            Uri.parse("package:" + getPackageName()))))
                    .setNegativeButton("Plus tard", null)
                    .show();
        } else if (checkSelfPermission(android.Manifest.permission.READ_EXTERNAL_STORAGE) != PackageManager.PERMISSION_GRANTED) {
            requestPermissions(new String[]{android.Manifest.permission.READ_EXTERNAL_STORAGE}, 0);
        }
    }

    void pickFolder() {
        Intent i = new Intent(Intent.ACTION_OPEN_DOCUMENT_TREE);
        i.putExtra(DocumentsContract.EXTRA_INITIAL_URI,
                DocumentsContract.buildDocumentUri("com.android.externalstorage.documents", "primary:Pictures"));
        startActivityForResult(i, PICK_FOLDER);
    }

    /** docId « primary:Pictures/x » -> /storage/emulated/0/Pictures/x ; « ABCD-1234:x » (carte SD) -> /storage/ABCD-1234/x. */
    static String treeToPath(String docId, String primaryRoot) {
        int c = docId.indexOf(':');
        String vol = docId.substring(0, c), rest = docId.substring(c + 1);
        String base = "primary".equals(vol) ? primaryRoot : "/storage/" + vol;
        return rest.isEmpty() ? base : base + "/" + rest;
    }

    @Override
    protected void onActivityResult(int req, int res, Intent data) {
        super.onActivityResult(req, res, data);
        if (req != PICK_FOLDER || res != RESULT_OK || data == null || data.getData() == null) return;
        String docId = DocumentsContract.getTreeDocumentId(data.getData());
        if (!docId.contains(":")) return;  // fournisseur autre que le stockage local (cloud…) : pas de chemin réel
        String path = treeToPath(docId, Environment.getExternalStorageDirectory().getAbsolutePath());
        // champ « dossier » et ouverture, comme si l'utilisateur avait saisi le chemin dans la page
        web.evaluateJavascript("document.getElementById('folder').value=" + JSONObject.quote(path) + ";openFolder();", null);
    }

    @Override
    public void onWindowFocusChanged(boolean focus) {
        super.onWindowFocusChanged(focus);
        if (focus) getWindow().getDecorView().setSystemUiVisibility(View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY
                | View.SYSTEM_UI_FLAG_FULLSCREEN | View.SYSTEM_UI_FLAG_HIDE_NAVIGATION
                | View.SYSTEM_UI_FLAG_LAYOUT_FULLSCREEN | View.SYSTEM_UI_FLAG_LAYOUT_HIDE_NAVIGATION);
    }

    @Override
    public void onBackPressed() {
        if (web.canGoBack()) web.goBack();
        else super.onBackPressed();
    }
}
