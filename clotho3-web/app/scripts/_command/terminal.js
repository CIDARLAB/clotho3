angular.module('clotho.commandbar')
.controller('TerminalCtrl', function($scope, ClothoCommandHistory) {
  $scope.logEntries = ClothoCommandHistory.entries;
});