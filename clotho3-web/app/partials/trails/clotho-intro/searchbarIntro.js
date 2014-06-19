$clotho.extensions.controller('clothoIntro_searchbarIntroCtrl', function($scope, $focus, $timeout, $modal, Clotho) {

	/*
    var dialogOpts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'The Clotho Command Bar',
                message: 'You can run all of Clotho\'s tools through the Command Bar. The Command Bar is capable of creating and editing your data, running commands, or searching. Type it in, and hit enter or the Submit Button.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    var dialog2Opts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'Entering Commands',
                message: 'You can enter commands in the text input. Submit by hitting enter, or clicking the Submit button.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    var dialog3Opts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'Activity Log',
                message: 'You can check your command and output history by opening the Activity Log. A counter tracks unread messages.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    $focus.addBackdrop();
    $focus.bringToFront($('#clothoSearchbar'))
    .then(function(oldZ) {
        return $timeout(function() { return oldZ }, 500 );
    })
    .then(function(oldZ) {
        return $modal.open(dialogOpts).result.then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#clothoSearchbar'));
    /*
        return $focus.bringToFront($('#searchBarInput').parent())
    })
    .then(function(oldZ) {
        return $timeout(function() { return oldZ }, 100 );
    })
    .then(function(oldZ) {
        return $modal.open(dialog2Opts).result.then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#searchBarInput').parent());
    *//*
        return $focus.bringToFront($('#searchbar_logbutton'))
    })
    .then(function(oldZ) {
        return $timeout(function() { return oldZ }, 100 );
    })
    .then(function(oldZ) {
        return $modal.open(dialog3Opts).result.then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#searchbar_logbutton'));
    })
    .then(function() {
        return $timeout(function() {$focus.removeBackdrop()}, 200);
    });
*/


    $scope.defineVariable = function() {
        var str = 'var sequence = "acgatcgatcgat"';
        $focus.typeOutSearch(str, true);
    };

    $scope.getObject = function() {
        var str = 'clotho.get(\'GFPmut3\')';
        $focus.typeOutSearch(str, true);
    };

    $scope.runQuery = function () {
        var str = 'clotho.query({"name":"GFPmut3"})';
        $focus.typeOutSearch(str, true);
    };

    $scope.runDNA = function() {
        var str = 'DNA.complement("acgttgca")';
        $focus.typeOutSearch(str, true);
    };

    //hidden:

    $scope.queryClothoTeam = function () {
        var str = 'clotho.query({schema : "ClothoTeam"})';
        $focus.typeOutSearch(str, true);
    };

    $scope.runFunction = function() {
        var str = 'clotho.run("uppercase", ["i will be ALL caps"])';
        $focus.typeOutSearch(str, true);
    };

});