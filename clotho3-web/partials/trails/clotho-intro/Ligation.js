'use strict';

Application.Extensions.controller('clothoIntro_LigationCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'DNA', 'Digest', 'PCR', function($scope, $focus, $timeout, $dialog, DNA, Digest, PCR) {
    $scope.DNA = DNA;
    $scope.Digest = Digest;
    $scope.PCR = PCR;

    $scope.ligate_demoSets = [
        ['aaaaaaaaaaA^CATG_', '^CATG_Tttggttggttgg'],
        ['aaaaaaaaaaA^CATG_', 'ccaaccaaccaaA^CATG_'],
        ['^CATG_Ttttttttttt', '^CATG_Tttggttggttgg'],
        ['^CATG_Ttttttttttt', 'ccaaccaaccaaA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'ggggggA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'gtcatcgatcagt_GTAC^'],
        ['aaaaaaaaaaA^CATG_T', 'A^CATG_Tacgatagcattaagcgt'],
        ['aaaaaaaaaaaaaa|', '|ttggttggttgg']
    ];

    $scope.setLigate = function(setInd) {
        $scope.fragments = $scope.ligate_demoSets[setInd];
    };

    $scope.setLigate(0);

    $scope.$watch(function() {
        return $scope.fragments[0] + $scope.fragments[1]
    }, function (newval, oldval) {
        console.log(newval);
        $scope.ligated = PCR.ligate($scope.fragments);
    });


    $scope.showMeHow = function() {

        //todo - move to funciton that takes array to type out commands

        $dialog.messageBox('Defining Variables', 'This time we\'ll define the variables as part of the process. First we need to define our two fragments, <code>frag1</code> and <code>frag2</code>. Then we\'ll join them in an array <code>fragments</code>, which we pass to the function <code>ligate</code>. <b>Remember we must pass our arguments in an array!</b>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
            .then(function() {
                var str = 'var frag1 = "aaaaaaaaaaA^CATG_";';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                var str = 'var frag2 = "^CATG_Tttggttggttgg";';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                var str = 'var fragments = [frag1, frag2];';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                return $focus.typeOutSearch('clotho.run("ligate", [fragments])');
            })
            .then(function() {
                $('#searchBarInput').focus()
            });
    };

}]);