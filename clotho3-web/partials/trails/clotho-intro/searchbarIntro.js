Application.Extensions.controller('clothoIntro_searchbarIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    var dialogOpts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: '/interface/templates/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'The Clotho Search Bar',
                message: 'All of Clotho\'s tools rely on lower-level functions, which you can run through the search bar if you need more control (or need to do things a little faster). Clotho\'s Search Bar is a powerful tool, able to do far more than search, and useable anywhere on the site. For example, you could create or edit your data, run commands, or, yes, even search.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    var dialog2Opts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: '/interface/templates/dialogMessagebox.html',
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
        templateUrl: '/interface/templates/dialogMessagebox.html',
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
        return $dialog.dialog(dialogOpts).open().then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#clothoSearchbar'));
        return $focus.bringToFront($('#searchBarInput').parent())
    })
    .then(function(oldZ) {
        return $timeout(function() { return oldZ }, 100 );
    })
    .then(function(oldZ) {
        return $dialog.dialog(dialog2Opts).open().then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#searchBarInput').parent());
        return $focus.bringToFront($('#searchbar_logbutton'))
    })
    .then(function(oldZ) {
        return $timeout(function() { return oldZ }, 100 );
    })
    .then(function(oldZ) {
        return $dialog.dialog(dialog3Opts).open().then(function() {return oldZ})
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#searchbar_logbutton'));
    })
    .then(function() {
        return $timeout(function() {$focus.removeBackdrop()}, 200);
    })



}]);