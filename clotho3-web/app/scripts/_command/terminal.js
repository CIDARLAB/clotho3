angular.module('clotho.commandbar').controller('TerminalCtrl', function($scope, CommandBar) {
    $scope.log = CommandBar.log;
});