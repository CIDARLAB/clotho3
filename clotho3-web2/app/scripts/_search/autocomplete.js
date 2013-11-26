angular.module('clotho.search').directive('commandBarAutocomplete', function(Clotho, CommandBar, $location, $window) {

	var commandTemplateUrl = 'views/_search/detail-command.html',
			authorTemplateUrl = 'views/_search/detail-author.html'

	return {
		restrict: 'A',
		templateUrl: 'views/_search/autocomplete.html'
	}


});