'use strict';

Application.Extensions.controller('clothoIntro_scriptingIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', function($scope, $focus, $timeout, $dialog) {

    $scope.showMeHow = function() {
        var searchInput = $('#searchBarInput'),
            searchSubmit = $('#searchBarSubmit'),
            oldinputZ = searchInput.css("z-index"),
            oldsubmitZ = searchSubmit.css("z-index"),
            maxZ = $focus.maxZ;


        $dialog.messageBox('Entering Commands', 'The search bar executes commands in JavaScript. You can run Clotho API commands, and functions within Clotho. For example, to reverse complement a sequence, you would call the function <code>revcomp</code> and pass your sequence in the array of arguments like so: <code>["acagtgcca"]</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function() {
            return $focus.typeOutSearch("clotho.run('revcomp', ['acagtgcca'])")
        })
        .then(function() {
            return $dialog.messageBox('Submit', 'You can submit your command by either hitting the enter key, or pressing the submit button', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        });

        /*//add backdrop, highlight input
        $focus.addBackdrop(maxZ+1)
        .then(function() {
            return $dialog.messageBox('Entering Commands', 'The search bar executes commands in JavaScript. You can run Clotho API commands, and functions within Clotho. For example, to reverse complement a sequence, you would call the function <code>revcomp</code> and pass your sequence in the array of arguments like so: <code>["acagtgcca"]</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        })
        .then(function() {
            return $focus.setZ(maxZ + 2, searchInput)
        })
        .then(function() {
            return $focus.typeOut(searchInput,
                "clotho.run('revcomp', ['acagtgcca'])", 'display.query')
        })
        .then(function() {
            return $timeout(function() {
                //fade out search
                $focus.setZ(oldinputZ, searchInput);
            }, 500)
        })
        .then(function() {
            return $dialog.messageBox('Submit', 'You can submit your command by either hitting the enter key, or pressing the submit button', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        })
        .then(function() {
            return $timeout(function() {
                $focus.removeBackdrop();
                $focus.setZ(oldsubmitZ, searchSubmit);
            })
        });
        .then(function() {
            return $dialog.messageBox('Submit', 'You can submit your command by either hitting the enter key, or pressing the submit button', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        })
        .then(function() {
            //focus submit
            return $focus.setZ(maxZ + 2, searchSubmit);
        })
        .then(function() {
            return $timeout(function() {
                return $dialog.messageBox('Getting your Result', 'The server will then execute your command. The server will communicate with you via the Activity Log, to the right of the search bar, and snippets of your most recent messages will appear for a few seconds. ', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
            }, 1500)
        })
        .then(function() {
            return $timeout(function() {
                //submit
                searchSubmit.click();

                $focus.removeBackdrop();
                $focus.setZ(oldsubmitZ, searchSubmit);
            })
        })
        .then(function() {
            return $timeout(function() {
                return $dialog.messageBox('Enter Your Own Command', 'Try running your own function! For example, you can transcribe a sequence of DNA by typing <code>clotho.run("transcribe", ["agctagctagcta"]</code>, or coming up with your own sequence.', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
            }, 1500)
        })
        .then(function() {
            searchInput.val('').focus();
        });*/
    };


    $scope.showHelpTips = function() {
        $focus.elementPopover('#searchBarInput', "whats up");
    };


}]);