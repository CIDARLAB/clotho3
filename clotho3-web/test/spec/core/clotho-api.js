//todo - need a way to better pass in options (e.g. socket) to provide a mock socket to test

describe('Clotho API', function() {
	var factory,
		socket,
		mockSocket,
		scope,
		spy;

	beforeEach(function() {

		module('clotho.core');

		inject(function(Clotho, $rootScope) {
			factory = Clotho;
			spy = jasmine.createSpy('callbackSpy');
			scope = $rootScope.$new();
		});

	});

	describe('#get', function () {

		it('should have a get function', function() {
			expect(angular.isFunction(factory.get)).toBe(true);
		});

	});

	describe('#set', function () {

		it('should have a set function', function() {
			expect(angular.isFunction(factory.set)).toBe(true);
		});

	});

});