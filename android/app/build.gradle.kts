import java.util.Properties

plugins {
    id("com.android.application")
    id("com.chaquo.python")
}

// Sources du dépôt (racine) recopiées à chaque build, sans les dupliquer dans git : modules Python sans les tests, et la page.
val repoRoot = rootProject.projectDir.parentFile
val gen = layout.buildDirectory.dir("petscan")
// Python 3.10 du poste (pip, compilation) : `buildPython=` dans local.properties, sinon Chaquopy cherche `py -3.10` / python3.10
val hostPython: String? = Properties().apply {
    rootProject.file("local.properties").takeIf { it.exists() }?.reader()?.use { r -> load(r) }
}.getProperty("buildPython")
val copyPython by tasks.registering(Sync::class) {
    from(repoRoot) { include("*.py"); exclude("test_*.py") }
    into(gen.map { it.dir("python") })
}
val copyPage by tasks.registering(Sync::class) {
    from(repoRoot) { include("index.html") }
    into(gen.map { it.dir("assets") })
}

// la copie des sources doit précéder tout le reste (sinon l'APK part sans le code Python ni la page)
tasks.named("preBuild") { dependsOn(copyPython, copyPage) }
tasks.configureEach { if (name != "copyPython" && name != "copyPage") mustRunAfter(copyPython, copyPage) }

android {
    namespace = "fr.petscan"
    compileSdk = 36
    buildToolsVersion = "36.0.0"
    defaultConfig {
        applicationId = "fr.petscan"
        minSdk = 26  // Chaquopy : 24 ; 26 pour le dossier de départ du sélecteur (EXTRA_INITIAL_URI)
        targetSdk = 34
        versionCode = 1
        versionName = "0.1.0"
        ndk { abiFilters += listOf("arm64-v8a", "x86_64") }  // tablettes + émulateur
    }
    sourceSets["main"].assets.srcDir(copyPage)
    compileOptions {
        sourceCompatibility = JavaVersion.VERSION_17
        targetCompatibility = JavaVersion.VERSION_17
    }
}

chaquopy {
    defaultConfig {
        version = "3.10"  // seule version pour laquelle le dépôt Chaquopy fournit scipy (1.8.1)
        hostPython?.let { buildPython(it) }
        pip {
            install("numpy")
            install("scipy")
            install("pydicom")
            install("matplotlib")
            install("pillow")
        }
    }
    sourceSets { getByName("main") { srcDir(copyPython) } }
}

