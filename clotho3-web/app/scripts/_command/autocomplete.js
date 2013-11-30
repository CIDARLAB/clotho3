//todo

angular.module('clotho.commandbar').directive('commandBarAutocomplete', function(Clotho, CommandBar, $location, $window) {

	var commandTemplateUrl = 'views/_command/detail-command.html',
			authorTemplateUrl = 'views/_command/detail-author.html'

	return {
		restrict: 'A',
		templateUrl: 'views/_command/autocomplete.html'
	}


});