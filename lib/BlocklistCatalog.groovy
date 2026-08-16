//
// Points to FooDMe2's built-in taxonomy blocklist.
//
// The revision is pinned so the same BarBeQuE version always uses the same
// blocklist, even when FooDMe2's main branch changes.
//
class BlocklistCatalog {

    static final String REVISION = 'ad8c8f854cd9a33d4edf995d42a8915eb73f71e8'
    static final String BASE_URL = "https://raw.githubusercontent.com/bio-raum/FooDMe2/${REVISION}"

    static String url() {
        return "${BASE_URL}/assets/blocklist.txt"
    }
}
