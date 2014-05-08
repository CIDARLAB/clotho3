angular.module('clotho.interface')
.run(function (hotkeys, $location, CommandBar) {
	hotkeys.add('f', 'Focus Command Bar', function (event) {
		event.preventDefault();
		CommandBar.focusInput();
	});
	hotkeys.add('a', 'Show Activity Log', function (event) {
		event.preventDefault();
		CommandBar.toggleActivityLog();
	});
	hotkeys.add('g h', 'Go to Homepage', function () {
		$location.path('/')
	});
	hotkeys.add('g b', 'Go to Browser', function () {
		$location.path('/browser')
	});
	hotkeys.add('g e', 'Go to Editor', function () {
		$location.path('/editor')
	});
	hotkeys.add('g t', 'Go to Trails', function () {
		$location.path('/trails')
	});
});