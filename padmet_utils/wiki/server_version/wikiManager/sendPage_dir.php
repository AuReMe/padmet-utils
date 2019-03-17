<?php

include 'globals.php';

$api_url = $argv[2];
$username = $argv[3];
$password = $argv[4];

function listFolderFiles($dir=null)

{
    $rtn = array();
    if (is_dir($dir)) {
        #echo '<ol>';
        foreach (new DirectoryIterator($dir) as $fileInfo) {
            if (!$fileInfo->isDot()) {
                if ($fileInfo->isDir()) {
                    $rtn = array_merge($rtn, listFolderFiles($fileInfo->getPathname()));
                } else {
                    $rtn[] = $fileInfo->getPathname();
                }
            }
        }
    } else {
        echo "No such directory '$dir'";
    }
    return $rtn;
}



try {
	$wiki = new Wikimate($api_url);
	if ($wiki -> login($username, $password))
		;
	else {
		$error = $wiki -> getError();
		echo "Wikimate error on login : " . $error['login'] . "\n";
		echo "Have to wait " . $error['wait'] . " seconds\n";
	}
} catch ( Exception $e ) {
	echo "<b>Wikimate error2</b>: " . $e -> getMessage() . "\n";
}

$files = listFolderFiles($argv[1]);
foreach($files as $file) {
	$filename = str_replace ("__47__", "/", basename($file,".txt"));
	echo $filename, "\n";
	$textpage = file_get_contents($file);
	// create a new page object
	$page = $wiki->getPage($filename);
	// check if the page exists or not and output the page
	$page->setText($textpage);
}


?>

