angular.module('clotho.search').directive('logEntries', function() {

	return {
		restrict: 'A',
		templateUrl: '../views/_search/logEntries.html',
		scope: {
			entries: '=logEntries'
		}
	}
});